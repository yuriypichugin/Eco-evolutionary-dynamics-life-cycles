#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 13:31:53 2020

@author: pichugin

Updated and cleaned collection of useful functions for homogeneous life cycle models 

!!! Is not backward compartible with v1 !!!

v2.0.1 added 'bumped_rand_exp' model
v2.1 added the second return to the InvasionRate() - condition number of the invasion matrix
	added 'random_bump_exponential' model to the BirthInit()
"""

"""
The life cycle data structure.

The life cycle is encoded as a list of pathways.
Each pathway is represented by a tuple: (Fragmented_Size, [Fragments_Sizes], Fragmentation_Rate, Trigger_Flag)

(int) Fragmented_Size - the size of the fragmented group
(int list) Fragments_Sizes - the list of size of fragments
(float) Fragmentation_Rate - the rate of fragmentation (for free pathways), or the probability of fragmentation (for triggered pathways)
(int) Trigger_Flag - the flag indicating the motive of the fragmentation

TriggerFlag==0 triggered pathways
TriggerFlag==1 free pathways (can be used by ProjectionMatrix(LifeCycle, MMSizes, b, d))

The sum of Fragments_Sizes should not exceed Fragmented_Size.
For a fragmentation without loss, this sum is equal to Fragmented_Size.

The sum of Fragmentation_Rate values for all triggered pathways with the same Fragmented_Size should not exceed one.
The sum of Fragmentation_Rate values for all triggered pathways with the maximal Fragmented_Size should be equal to one.
"""

import numpy as np
from IntegerPartitions import revlex_partitions
from copy import deepcopy as deepcopy
from scipy.integrate import solve_ivp
from numpy.linalg import cond

"""
SECTION 1:
Initialization of environmental conditions
"""

def BirthInit(model, NMax, Params):
    """
    Initialization of the birth rates vector
    
    Arguments:
    (string) model - the model used
    (int) NMax  - the size of vector (the maximal size of group existing in a population)
    (list of floats) Params - parameters of the birth rates model
    
    Return:
    (list of floats) b - the vector of birth rates (per group size!)
    """
    
    b=[0]*NMax
    
    if model=='direct':
        """rates are equal to Params"""
        for i in range(NMax):
            b[i] = Params[i]
            
    elif model=='power':
        """ Power model. """
        M = Params[0]
        S = Params[1]
        for i in range(NMax):
            i_eff = i/(NMax-1.0)
            b[i] = 1 + M*pow(i_eff,S)
            
    elif model=='rand_uniform':
        """ Random uniform model. """
        for i in range(NMax):
            b[i] = np.random.uniform(0.0, 1.0)
            
    elif model=='rand_exp':
        """ Random exponential model """
        Lambda = Params[0]
        for i in range(NMax):
            b[i] = np.random.exponential(Lambda)
            
    elif model=='rand_norm':
        """ Random gauss model """
        Mean = Params[0]
        Std = Params[1]
        for i in range(NMax):
            b[i] = np.random.normal(Mean, Std)
            
    elif model == 'rand_sorted_up':
        """ monotonically increasing random sequence """
        true_bi = np.random.uniform(0.0, 1.0, NMax)
        true_bi_sorted = np.sort(true_bi)
        for i in range(NMax):
            b[i] = true_bi_sorted[i]
            
    elif model == 'rand_sorted_down':
        """ monotonically decreasing random sequence """
        true_bi = np.random.uniform(0.0, 1.0, NMax)
        true_bi_sorted = np.sort(true_bi)[::-1]
        for i in range(NMax):
            b[i] = true_bi_sorted[i]
            
    elif model == 'rand_sorted_unimodal':
        """ unimodal random sequence """
        LocMax = Params[0]
        true_bi = np.random.uniform(0.0, 1.0, NMax)
        true_bi_sorted = np.sort(true_bi)[::-1]
        LeftSlot = LocMax - 1 
        RightSlot = LocMax + 1 
        unimodal_true_bi = np.zeros(NMax)
        unimodal_true_bi[LocMax] = true_bi_sorted[0]
        b[LocMax] = unimodal_true_bi[LocMax]
        for i in range(1, NMax):
            if (LeftSlot >= 0) and (RightSlot < NMax):
                coin = np.random.randint(2)
                if coin == 0:
                    position = LeftSlot
                    LeftSlot -= 1
                else:
                    position = RightSlot
                    RightSlot += 1
            elif (LeftSlot < 0) and (RightSlot < NMax):
                position = RightSlot
                RightSlot += 1
            elif (LeftSlot >= 0) and (RightSlot >= NMax):
                position = LeftSlot
                LeftSlot -= 1
            else:
                position = 0
                print(" profile generation went wrong ! ")
            unimodal_true_bi[position] = true_bi_sorted[i]
            b[position] = unimodal_true_bi[position]

    elif model=='bumped_rand_exp':
        """ Random exponential model plus bump """
        Lambda = Params[0]
        Bump = Params[1]
        for i in range(NMax):
            b[i] = np.random.exponential(Lambda) + Bump

            
    else:
        """ default model, no effect of the group size on the growth of cells """
        for i in range(NMax):
            b[i]=1.0
    return b 

def DeathInit(model, NMax, Params):
    """
    Initialization of the death rates vector
    
    Arguments:
    
    (string) model - the number of model used
    (int) NMax  - the size of vector (the maximal size of group existing in a population)
    (list of floats) Params - parameters of the fitness landscape
    
    Return:
    (list of floats) d - the vector of death rates
    """
    d=[0]*NMax
    
    if model == 'direct':
        """Small sizes scans model. Death is determined by params """
        for i in range(NMax):
            d[i] = Params[i]

    elif model == 'power':
        """ Power model. """
        M = Params[0]
        S = Params[1]
        for i in range(NMax):
            i_eff = i/(NMax-1.0)
            d[i] = -M*pow(i_eff,S)
            
    elif model == 'rand_uniform':
        """ Random uniform model. """
        for i in range(NMax):
            d[i] = np.random.uniform(0.0, 1.0)
            
    elif model == 'rand_exp':
        """ Random exponential model """
        Lambda = Params[0]
        for i in range(NMax):
            d[i] = np.random.exponential(Lambda)
            
    elif model == 'rand_norm':
        """ Random gauss model """
        Mean = Params[0]
        Std = Params[1]
        for i in range(NMax):
            d[i] = np.random.normal(Mean, Std)
            
    elif model == 'rand_sorted_up':
        """ monotonically increasing random sequence """
        di = np.random.uniform(0.0, 1.0, NMax)
        di_sorted = np.sort(di)
        for i in range(NMax):
            d[i] = di_sorted[i]
            
    elif model == 'rand_sorted_down':
        """ monotonically decreasing random sequence """
        di = np.random.uniform(0.0, 1.0, NMax)
        di_sorted = np.sort(di)[::-1]
        for i in range(NMax):
            d[i] = di_sorted[i]
            
    elif model == 'rand_sorted_unimodal':
        """ unimodal random sequence """
        LocMax = Params[0]
        di = np.random.uniform(0.0, 1.0, NMax)
        di_sorted = np.sort(di)
        LeftSlot = LocMax - 1 
        RightSlot = LocMax + 1 
        unimodal_di = np.zeros(NMax)
        unimodal_di[LocMax] = di_sorted[0]
        d[LocMax] = unimodal_di[LocMax]
        for i in range(1, NMax):
#            print([LeftSlot, RightSlot])
            if (LeftSlot >= 0) and (RightSlot < NMax):
                coin = np.random.randint(2)
                if coin == 0:
                    position = LeftSlot
                    LeftSlot -= 1
                else:
                    position = RightSlot
                    RightSlot += 1
            elif (LeftSlot < 0) and (RightSlot < NMax):
                position = RightSlot
                RightSlot += 1
            elif (LeftSlot >= 0) and (RightSlot >= NMax):
                position = LeftSlot
                LeftSlot -= 1
            else:
                position = 0
                print(" profile generation went wrong ! ")
#            print(position)
            unimodal_di[position] = di_sorted[i]
            d[position] = unimodal_di[position]
            
    else:
        """ default model, no effect of the group size on the death of cells """
        for i in range(NMax):
            d[i] = 0
    return d


def InteractionInit(model, NMax, Params):
    """
    Initialization of the interaction rates matrix
    
    Arguments:
    
    (string) model - the number of model used
    (int) NMax  - the size of vector (the maximal size of group existing in a population)
    (list of floats) Params - parameters of the fitness landscape
    
    Return:
    (numpy array) K - the matrix of rates
    """
    K = np.zeros((NMax, NMax))
    if model == 'const':
        C = Params[0]
        K = C * np.ones((NMax, NMax))
        
    elif model == 'singleton':
        Value = Params[0]
        Killer = Params[1]-1
        Victim = Params[2]-1
        K[Victim, Killer] = Value
    
    elif model == 'kernel_killer':
        for i in np.arange(NMax):
            K[:, i] = Params[i]

    elif model == 'kernel_victim':
        for i in np.arange(NMax):
            K[i, :] = Params[i]

    elif model == 'rand_uniform':
        K = np.random.uniform(0.0, 1.0, (NMax, NMax))
    
    elif model == 'rand_exp':
        Lambda = Params[0]
        K = np.random.exponential(Lambda, (NMax, NMax))

    elif model == 'bumped_rand_exp':
        Lambda = Params[0]
        Bump = Params[1]
        K = np.random.exponential(Lambda, (NMax, NMax)) + Bump
     
    elif model == 'direct':
        for i in np.arange(NMax):
            for j in np.arange(NMax):
                K[i,j] = Params[i][j]
    else:
        K = np.ones((NMax, NMax))
    return K













"""
SECTION 2:
Generation of life cycle structures (see the top of the file for the structure definition)
"""

def PartList(N):
    """
    Calculation of the partition list in a fragmentation without losses.
    The sum of offspring sizes is equal to the size of the fragmented group
    
    Arguments:
    (int) N - the maximal size of group presented in population.
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for Size in range(2, N+2): 
        #create generator of partitions of Size (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(Size)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            LifeCycle.append([(Size, PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    LC_ID=[]
    for counter in range(len(LifeCycle)):
        LC_ID.append(LifeCycle[counter])

    return LC_ID


def PartList_Loss(N, Loss):
    """
    Calculation of the partition list in a fragmentation with given losses.
    This is intended to replace PartList_WL().
    The sum of offspring sizes is equal to the size of the fragmented group minus Loss
    
    Arguments:
    (int) N - the maximal size of group presented in population.
    (int) Loss - the number of cells lost in a fragmentation
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for SizeOff in range(2, N+2-Loss): 
        #create generator of partitions of SizeOff (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(SizeOff)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            LifeCycle.append([(SizeOff+Loss, PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    LC_ID=[]
    for counter in range(len(LifeCycle)):
        LC_ID.append(LifeCycle[counter])

    return LC_ID

def PartList_ChainLoss(N, Loss):
    """
    Calculation of the partition list in a fragmentation with given proportional losses.
    This case represent the fragmentation of linear filaments
    The sum of offspring sizes is equal to the size of the fragmented group minus Loss*(#parts-1)
    
    Arguments:
    (int) N - the maximal size of group presented in population.
    (int) Loss - the number of cells lost at each break of the chain
    NOTE: the partition goes up to sizes N+1, since the fragmentation follows the growth of group of size N into size N+1
    
    Return:
    (list of lists) LC_ID - the list of partitions. Each partition contain its ID number and the LifeCycle tuple that describes the fragmentation pattern
    """
    #create list of life cycles
    LifeCycle=[]
    for SizeOff in range(2, N+2): 
        #create generator of partitions of SizeOff (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(SizeOff)
        #convert generator into a list of partititons
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        #Transform list of partitions into list of life cycles
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            Num_of_parts=len(PartList[Reaction])            
            LifeCycle.append([(SizeOff+Loss*(Num_of_parts - 1), PartList[Reaction], 1, 0)])
    
    #attach an ID number to each life cycle
    pre_LC_ID=[]
    for counter in range(len(LifeCycle)):
        pre_LC_ID.append(LifeCycle[counter])

    LC_ID=[]
    for ID in pre_LC_ID:
        if ID[1][0][0]<(N+2):
            LC_ID.append(ID)
    
    return LC_ID

def Generate_Random_Stochastic_LC(Nmax):
    """
    Generates life cycle as a list of tuples
    
    Arguments:
    Nmax - the maximal size of groups EXISTING in the population
    
    Return:
    LifeCycle - life cycle structure
    """
    TrigProbMagnitude=np.random.uniform(0.0, 1.0, Nmax+1)     #the sum of probabilities of all rates of triggered fragmentation occurring with groups of a each size
    
    LifeCycle=[]
    #set up triggered reactions
    for Size in range(2, Nmax+2): 
        #create generator of partitions of Size (i.e. list of fragments sizes)
        PartitionGenerator=revlex_partitions(Size)
        #convert generator into a list 
        PartList=[]
        for Part in PartitionGenerator:
            L=len(Part)
            E=[0]*L
            for j in range(L):
                E[j]=Part[j]
            PartList.append(E)
        
        #compute fragmentation probabilities
        Probs=np.random.exponential(1, len(PartList))
        if Size < Nmax+1:
            Norm=TrigProbMagnitude[Size]/sum(Probs[1:])
        else:                                           #the combined fragmentation probability at the maximal size is equal to one
            Norm=1/sum(Probs[1:])
        Probs=Probs*Norm
        
        #Construct biological reactions data
        for Reaction in range(1,len(PartList)): #the first partition is trivial and therefore does not correspond to a reaction
            RData=(Size, PartList[Reaction], Probs[Reaction], 0)        
            LifeCycle.append(RData)
    
    return LifeCycle


def Generate_Mixed_LC(PureLC_1, PureLC_2, Fraction_1):
    """
    Generates a mixed life cycle with two partitions
    
    Arguments:
        (Life cycle structure) PureLC_1 - the first pure life cycle in the mix
        (Life cycle structure) PureLC_2 - the second pure life cycle in the mix
        (float) Fraction_1 - the probability for the PureLC_1 to be executed
        Both life cycles must have only one partition each!
                
    Return:
        (Life cycle structure) OutputLC - mixed life cycle
    """
    
    if len(PureLC_1)>1:
        print("The first component is not a pure life cycle!")
    if len(PureLC_2)>1:
        print("The second component is not a pure life cycle!")
    
    MaturitySize_1 = MinMaxSizes(PureLC_1)[1]
    MaturitySize_2 = MinMaxSizes(PureLC_2)[1]
    
    MixedLC = [list(PureLC_1[0]), list(PureLC_2[0])]
    
    if MaturitySize_1 > MaturitySize_2:
        MixedLC[0][2] = 1.0
        MixedLC[1][2] = 1.0 - Fraction_1
        
    if MaturitySize_1 == MaturitySize_2:
        MixedLC[0][2] = Fraction_1
        MixedLC[1][2] = 1.0 - Fraction_1
        
    if MaturitySize_1 < MaturitySize_2:
        MixedLC[0][2] = Fraction_1
        MixedLC[1][2] = 1.0
  
    OutputLC = [tuple(i) for i in MixedLC]

    return OutputLC











"""
SECTION 3:
Construction of the projection matrix and finding exponential growth rates
"""

def MinMaxSizes(LifeCycle):
    """
    calculates the minimal and maximal sizes of groups 
    that appears in the life cycle.
    This is later used for the projection matrix construction
       
    Arguments:
    (list of tuples) LifeCycle - life cycle data, see description of the format in the separate comment
       
    Return:
    The list of two integers: the first one is the minimal size of groups, the second one is the maximal size of groups in the course of a given life cycle
    """
    MaxSize=0
    for pathway in range(len(LifeCycle)):
        Fragmenting_Group_Size = LifeCycle[pathway][0]
        TriggerFlag = LifeCycle[pathway][3]
        if MaxSize < (Fragmenting_Group_Size - 1 + TriggerFlag):
            MaxSize = (Fragmenting_Group_Size - 1 + TriggerFlag)
    
    MinSize=MaxSize
    for pathway in range(len(LifeCycle)):
        Offsprings=LifeCycle[pathway][1]
        MinOffspring=min(Offsprings)
        if MinOffspring<MinSize:
            MinSize=MinOffspring
        
    return [MinSize, MaxSize]


def ProjectionMatrix(LifeCycle, b, d, MatrixSize):
    """
    calculates the projection matrix of a given life cycle
    
    Arguments:
    
    (list of tuples) LifeCycle - life cycle data, see description of the format in the separate comment
    (list of floats) b - the list of the growth rates
    (list of floats) d - the list of the death rates
    (int) MatrixSize - the size of the resulting matrix
    
    Return:
    (matrix) A - projection matrix corresponding to a given life cycle in a given fitness landscape
    """
    
#    MMSizes=MinMaxSizes(LifeCycle)
#    MatrixSize=MMSizes[1]
    
    A=np.zeros((MatrixSize, MatrixSize))
    A=np.asmatrix(A)
    
    #set up the fitness landscape
    for i in range(MatrixSize):
        A[i,i] = A[i,i] - (i+1)*b[i]
        A[i,i] = A[i,i] - d[i]
        if i < (MatrixSize-1):
            A[i+1,i] = A[i+1,i] + (i+1)*b[i]
    
    #include the impact of fragmentation pathways
    for pathway in range(len(LifeCycle)):
    
        Pathway=LifeCycle[pathway]
    
        if Pathway[3]==0:        #triggered pathways
            Rate=Pathway[2] * (Pathway[0] - 1) * b[Pathway[0] - 2]
            Column=Pathway[0] - 2 
        elif Pathway[3]==1:      #free pathways
            Rate=Pathway[2]
            Column=Pathway[0] - 1 
        
        for i in range(len(Pathway[1])):
            Row=Pathway[1][i] - 1 
            A[Row, Column] = A[Row, Column] + Rate
        
        if Pathway[3]==0:
            if Column < (MatrixSize-1):
                A[Column+1,Column]=A[Column+1,Column] - Rate
        elif Pathway[3]==1:
            A[Column,Column]=A[Column,Column] - Rate
        
    return A



def ExponentialGrowthRate(LifeCycle, b, d):
    """
    calculates the proliferation rate provided by a life cycle under given fitness landscape
    
    Arguments:
    
    (list of tuples) LifeCycle - life cycle data, see description of the format elsewhere
    (list of floats) b - the list of the growth rates
    (list of floats) d - the list of the death rates
    
    Return:
    (float) ProliferationRate - proliferation rate provided by the life cycle in a given fitness landscape
    """    
       
    #Find the minimal and the maximal size of groups emerging in the given life cycle
    MMSizes=MinMaxSizes(LifeCycle)
    
    #check 1: the length of b and d vectors should be at least MaxSize 
    if (len(b)<MMSizes[1]):
        print ("Error: vector b is too short for the life cycle")
        return float('NaN')
    if(len(d)<MMSizes[1]):
        print("Error: vector d is too short for the life cycle")
        return float('NaN')
        
        
    
    #Initialize the projection matrix
    A=ProjectionMatrix(LifeCycle, b, d, MMSizes[1])
    
    #Calculate the proliferation rate
    Eigenvalues=np.linalg.eig(A)[0]
    ProliferationRate=max(Eigenvalues.real)

    return ProliferationRate










"""
SECTION 4:
Eco-evo dynamics
"""

def InvasionRate(A, K, X):
    """
    calculates the invasion rate in a given population
    
    Arguments:
    
    (matrix)         A - the projection matrix
    (numpy array)    K - the matrix of interaction rates
    (numpy array)    X - the current population composition
    
    Return:
    (float) InvasionRate - the rate of invasion into current population
    """
    DiffMatrix = A - np.diag( np.dot(K, X) )
    Eigenvalues=np.linalg.eig(DiffMatrix)[0]
    InvasionRate=max(Eigenvalues.real)
    
    Cond_number = cond(DiffMatrix)
	
    return [InvasionRate, Cond_number]
    
    

def StationaryState(LifeCycle, b, d, K, EndMaxSize):
    """
    calculates the equilibrium demographic distribution of a single life cycle
    
    Arguments:
    
    (list of tuples) LifeCycle - life cycle data, see description of the format on top of the file
    (list of floats) b - the list of the growth rates
    (list of floats) d - the list of the death rates
    (numpy array)    K - the matrix of interaction rates
    (int)            EndMaxSize - the size of the final vector
    
    Return:
    (matrix) X_equilibrium - equilibrium distribution of group sizes
    """
    
    """ constants """
    T_record = 1000
    
    """ compute the projection matrix """
    MMSizes=MinMaxSizes(LifeCycle)
    MaxSize=MMSizes[1]
    A = ProjectionMatrix(LifeCycle, b, d, MaxSize)
    K = K[0:MaxSize, 0:MaxSize]

    """ set up the initial approximation of X0 """
    MaxSize = MinMaxSizes(LifeCycle)[1]
    X0 = np.ones(MaxSize)
    meanKX = np.mean(np.asarray( np.dot( K, X0 ) ))
    LeadEV = max(np.linalg.eig(A)[0].real)
    X0 = X0 * LeadEV/meanKX
    
    """ simulate the dynamics to get the equilibrium """
    Record = LifeCyclesCompetition([A], K, [X0], [T_record])
    X_equilibrium = Record[0,0,:]
    
#    Norm_dX_target = 1e-5
#    Scale = 0.01
#
#    DiffMatrix = A - np.diag( np.dot(K, X0) )
#    dX = np.asarray( np.dot( DiffMatrix, X0 ) )[0]
#    Norm_dX = np.sum(np.abs(dX))/len(dX)
#        
##    step = 0
#    while Norm_dX > Norm_dX_target:
# 
#        X0 +=  Scale * dX
#        DiffMatrix = A - np.diag( np.dot(K, X0) )
#        dX = np.asarray( np.dot( DiffMatrix, X0 ) )[0]
#        Norm_dX = np.sum(np.abs(dX))/len(dX)

#        step+=1
#    print(step)
    """ update the length of the result to the standard """
    if len(X_equilibrium)<EndMaxSize:
        X_equilibrium = np.append(X_equilibrium, np.zeros(EndMaxSize-len(X_equilibrium)))
    
    return X_equilibrium

def PopDynamics(t, X, A, K):
    """
    computes the derivative of the population state.
    used to define ODE.
    """
    DiffMatrix = A - np.diag( np.dot(K, X) )
    dX = np.asarray(np.dot( DiffMatrix, X ))[0]
    return dX

def LifeCyclesCompetition(ProjectionMatricesSet, InteractionsMatrix, InitialStates, RecordTimes):
    """
    Calculates the dynamics of competing populations with different life cycles
    
    Arguments:
        (list of 2d arrays) ProjectionMatricesSet        - the list of projection matrices of the same size
        (2d array) InteractionsMatrix                    - the matrix of interactions (of the same size)
        (list of 1d arrays) InitialStates                - the list of population compositions for each of populations
        (list) RecordTimes                               - the list of time steps at which the population state must be recorded
    
    return:
        (3d numpy array) PopulationRecord      - the structure in the format
                                                            PopulationRecord[i,:,:] is the record of i-th subpopulation dynamics
                                                            PopulationRecord[i,j,:] is the state of i-th subpopulation at the moment RecordTimes[j]
                                                            PopulationRecord[i,j,k] is the number of groups of size k+1 in i-th population at the moment RecordTimes[j]
    """
    
    """ check inputs size agreement """
    Num_projection = len(ProjectionMatricesSet)
    Num_states = len(InitialStates)
    if Num_projection == Num_states:
        CompetitorsNum = Num_projection
    else:
        print("Numbers of projection matrices and initial states are not the same")
        return 0
    
    Sizes_projection = list(map(len, ProjectionMatricesSet))
    Sizes_interaction = len(InteractionsMatrix)
    Sizes_states = list(map(len, InitialStates))
    
    ArraySize = np.max(Sizes_projection)
    if not(ArraySize == np.min(Sizes_projection)):
        print("Projection matrices have different sizes")
        return 0
    if not(np.max(Sizes_states) == np.min(Sizes_states)):
        print("Initial states have different lengths")
        return 0
    if not(ArraySize == Sizes_interaction):
        print("Interaction and projection matrices have different sizes")
        return 0
    if not(ArraySize == np.max(Sizes_states)):
        print("Initial states and projection matrices have different sizes")
        return 0
    
    """ construct combined projection matrix """
    Combined_projection = []
    ZeroMatrix = np.zeros((ArraySize, ArraySize))
    for i in np.arange(CompetitorsNum):
        OneRow = [ZeroMatrix]*CompetitorsNum
        OneRow[i]=ProjectionMatricesSet[i]
        Combined_projection.append(OneRow)
    Combined_projection = np.block(Combined_projection)
    
    """ construct combined interaction matrix """
    Combined_interaction = []
    for i in np.arange(CompetitorsNum):
        OneRow = [InteractionsMatrix]*CompetitorsNum
        Combined_interaction.append(OneRow)
    Combined_interaction = np.block(Combined_interaction)
    
    """ construct combined population """
    X_init = np.concatenate(InitialStates)
    
    """ perform integration """
#    Derivative = PopDynamics([], X_init, Combined_projection, Combined_interaction)
    DynSol = solve_ivp(PopDynamics, [0, np.max(RecordTimes)], X_init, t_eval = RecordTimes, args=(Combined_projection, Combined_interaction)).y
    
    """ format output """
    PopultionRecord = []
    for i in np.arange(CompetitorsNum):
        Competitor_record = np.transpose(DynSol[(ArraySize*i):(ArraySize*(i+1)),:])
        PopultionRecord.append(Competitor_record)
    PopultionRecord = np.asarray(PopultionRecord)
    
    return PopultionRecord