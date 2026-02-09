import numpy as np
import pandas as pd 
import random
import math 
from itertools import product
from itertools import combinations_with_replacement
import matplotlib.pyplot as plt
from tqdm import tqdm

#------------------------------------
#   generating the genotype dictinary
#------------------------------------
def  two_locus_dictionary_function(allele_number = 3):

    haplotype = [
        f'A{i+1} B{j+1}'
        for i, j in product(range(allele_number), repeat=2)
    ]

    counter = 0
    genotype_dic = {}
    
    for h1, h2 in combinations_with_replacement(haplotype, 2): 
        genotype_dic[f'{h1} / {h2}'] = counter
        genotype_dic[f'{h2} / {h1}'] = counter
        counter += 1 
    return genotype_dic 

def genotype_converter (genotype, genotype_dictionary):
    return genotype_dictionary.get(genotype, 1000)

#----------------------------
#   generating the population
#----------------------------
def genotype_generator (): 

    genotype_dic = two_locus_dictionary_function()
    index = np.random.randint(0, len(genotype_dic))
    genotype = [*genotype_dic.keys()][index]

    return genotype

def population_dataframe_generator(IP): 
    return pd.DataFrame({
        "Mating type": ["C / C"] * IP,
        "Genotype": [genotype_generator() for _ in range(IP)],
        "Parasite status": np.random.choice(
            [True, False],
            size = IP
            )
    })

#--------------------------
#   calculating frequencies
#--------------------------
def parasite_genotype_frequencies(herm_population, male_population, genotype_dic):


    population = pd.concat((herm_population, male_population), ignore_index=True)
    infected = population[population['Parasite status'] == True]

    if len(infected) == 0: 
        return [0.0] * len(set(genotype_dic.values()))
    
    alpha = [0.0] * len(set(genotype_dic.values()))
    for genotype in infected['Genotype']: 
        u = genotype_converter(genotype, genotype_dic)
        alpha[u] += 1/len(infected)

    return alpha

def mating_type_frequency(s_frequency, c_frequency, herm_population, male_population):

    allele_map = {
        'C / C' : (1.0, 0.0), 
        'C / S' : (0.5, 0.5),
        'S / C' : (0.5, 0.5),
        'S / S' : (0.0, 1.0)
    }

    population = pd.concat([male_population, herm_population], ignore_index=True)

    N = len(population)
    c_count = 0
    s_count = 0 
    
    for mt in population['Mating type']: 
        c, s = allele_map[mt]
        c_count += c
        s_count += s
    try: 
        s_frequency.append(s_count / N)
        c_frequency.append(c_count / N)
    except ZeroDivisionError: 
        if c_count == 0:
            s_frequency, c_frequency = (1, 0)
        else: 
            s_frequency, c_frequency = (0, 1)
    except exception as e: 
        return None

    return s_frequency, c_frequency  

#---------------------------------
#   death  and infection function. 
#---------------------------------
def death_function (population, individual_index):
    population = population.drop(individual_index).reset_index(drop = True)
    return population

def infection_function (population, individual_index):
    population.at[individual_index , 'Parasite status'] = False \
        if population.at[individual_index, 'Parasite status'] else True
    return population

#----------------------------------------
#   functions required for birth function
#----------------------------------------
def mating_male_index (male_population, infection_effect_on_male): 
    weights = []
    for status in male_population['Parasite status']:
        if status: 
            weights.append(1-infection_effect_on_male)
        else:
            weights.append(1.0)
    weights = np.array(weights)
    weights /= weights.sum()

    return np.random.choice(len(male_population), p=weights)

def recombination (gen, r):
    if isinstance(gen, str):
        gen = gen.split()

    if np.random.rand() < r:
        gen[1], gen[4] = gen[4], gen[1]

    return gen

def newborn_genotype_generator(herm_gen, male_gen, r):
    herm_gen = recombination(herm_gen, r)
    male_gen = recombination(male_gen, r)

    match np.random.randint(2):
        case 0:
            idx_sperm = (0, 1)
        case 1:
            idx_sperm = (3, 4)
    
    match np.random.randint(2): 
        case 0: 
            idx_ovule = (0, 1)
        case 1: 
            idx_ovule = (3, 4)

    sperm_genotype = " ".join(male_gen[idx_sperm[0] : idx_sperm[1]+1])
    ovule_genotype = " ".join(herm_gen[idx_ovule[0] : idx_ovule[1]+1])

    return f'{sperm_genotype} / {ovule_genotype}'

def newborn_gender_generator (r):
    return 'herm' if np.random.rand() < r else 'male'

def mating_type_generator (herm_population,
                           individual_index,
                           male_population,
                           index_of_mating_male):
    
    mt_ovule = herm_population.at[individual_index, 'Mating type'].split()
    mt_ovule = mt_ovule[2*np.random.randint(2)]

    mt_sperm = male_population.at[index_of_mating_male, 'Mating type'].split()
    mt_sperm = mt_sperm[2*np.random.randint(2)]

    return f'{mt_sperm} / {mt_ovule}'

#-----------------
#   birth function
#-----------------
def birth_function (herm_population, male_population, individual_index, r,
                    infection_effect_on_male = 0, ro = 0.5, rs = 0.9, o_r = 0.5):

    male_index = 0
    newborn_parasite_status = False
    herm_gen = herm_population.at[individual_index, 'Genotype']
    herm_mt = herm_population.at[individual_index , 'Mating type']

    if herm_mt == "C / C": 
        outcross = True
        gender_prob = ro
    elif herm_mt in ["C / S", "S / C"]:
        outcross = np.random.rand() < o_r
        gender_prob = ro if outcross else rs
    elif herm_mt == "S / S": 
        outcross = False
        gender_prob = rs

    if outcross: 
        male_index = mating_male_index(
            male_population,
            infection_effect_on_male
        )
        male_gen = male_population.at[male_index, 'Genotype']
        mt = mating_type_generator(herm_population,
                                   individual_index,
                                   male_population,
                                   male_index
                                )
    else:
        male_index = individual_index
        male_gen = herm_gen
        mt = mating_type_generator(
            herm_population,
            individual_index,
            herm_population,
            individual_index
            )

    newborn_genotype = newborn_genotype_generator(herm_gen, male_gen, r)
    newborn_gender = newborn_gender_generator(gender_prob)

    if newborn_gender == 'male':
        male_population.loc[len(male_population)] = [
            mt,
            newborn_genotype, 
            newborn_parasite_status
        ]

    if newborn_gender == 'herm':
        herm_population.loc[len(herm_population)] = [
            mt,
            newborn_genotype,
            newborn_parasite_status
        ]

    return herm_population, male_population


#-----------------------------------
#   Gillespie Algorithm requiremnets
#-----------------------------------
def tau_generator(a): 
    if a <= 0:
        return math.inf
    return np.random.exponential(1/a)

def next_event(ta_herm, ta_male): 
    min_tau_herm = np.min(ta_herm) if len(ta_herm) > 0 else math.inf
    min_tau_male = np.min(ta_male) if len(ta_male) > 0 else math.inf

    if min_tau_herm <= min_tau_male: 
        tau = min_tau_herm
        idx_event = np.argmin(ta_herm)
        sex = 'Herm'
        reaction_type = REACTION_MAP['herm'].get(idx_event % 3, "Unkown")
        individual_index = idx_event // 3
    else: 
        tau = min_tau_male
        idx_event = np.argmin(ta_male)
        sex = 'Male'
        reaction_type = REACTION_MAP['male'].get(idx_event % 2, "Unkown")
        individual_index = idx_event // 2
    return tau, sex, reaction_type, individual_index

#-----------------------------------------------------
#   running population simulation for t_max time unit
#-----------------------------------------------------
def simulate_population(**kwargs):

    herm_population = kwargs['hp'].copy(deep=True)
    male_population = kwargs['mp'].copy(deep=True)

    t = 0
    iteration = 0
    time_series = [t]
    herm_sizes = [len(herm_population)]
    male_sizes = [len(male_population)]
    reactions = []

    s_frequency = []
    c_frequency = []
    s_frequency , c_frequency = mating_type_frequency(s_frequency,
                                                       c_frequency,
                                                       herm_population,
                                                       male_population
                                                       )

    if kwargs['scheduled_events'] is None:
        kwargs['scheduled_events'] = {}

    alpha = kwargs['alpha']

    while t < kwargs['t_max']:
        herm_population['Genotype Index'] = herm_population['Genotype'].apply(
        lambda g: genotype_converter(g, kwargs['genotype_dictionary'])
        )
        male_population['Genotype Index'] = male_population['Genotype'].apply(
        lambda g: genotype_converter(g, kwargs['genotype_dictionary'])
        )
        if iteration in kwargs['scheduled_events']:
            herm_population.at[0, 'Mating type'] = 'S / S'
        total_pop = len(herm_population) + len (male_population)

        ta_herm = []
        for __,row in herm_population.iterrows():
            g_index = row['Genotype Index']
            if row['Parasite status']: 
                ta_herm.extend([
                    tau_generator(kwargs['d']),
                    tau_generator(kwargs['b0i'] - (kwargs['ai'] * total_pop)),
                    tau_generator(kwargs['beta'])
                ])
            else: 
                ta_herm.extend([
                    tau_generator(kwargs['d']),
                    tau_generator(kwargs['b0u'] - (kwargs['au'] * total_pop)),
                    tau_generator(alpha[g_index])
                ])
        ta_male = []
        for __,row in male_population.iterrows():
            g_index = row['Genotype Index']
            if row['Parasite status']: 
                ta_male.extend([
                    tau_generator(kwargs['d']),
                    tau_generator(kwargs['beta'])
                ])
            else: 
                ta_male.extend([
                    tau_generator(kwargs['d']),
                    tau_generator(alpha[g_index])
                ])

        tau, sex, reaction_type, individual_index = next_event(ta_herm, ta_male)
        herm_population = herm_population.drop(columns = ['Genotype Index'])
        male_population = male_population.drop(columns = ['Genotype Index'])
        if sex == 'Herm': 
            match reaction_type: 
                case 'Birth':
                    herm_population, male_population = birth_function(herm_population,
                                                                    male_population,
                                                                    individual_index,
                                                                    kwargs['r'],
                                                                    kwargs['infection_effect_on_male'],
                                                                    kwargs['ro'],
                                                                    kwargs['rs'],
                                                                    kwargs['o_r'])

                case 'Death': 
                    herm_population = death_function(herm_population, individual_index)
                case 'Infection / Recovery': 
                    herm_population = infection_function(herm_population, individual_index)
        elif sex == 'Male': 
            match reaction_type: 
                case 'Death': 
                    male_population = death_function(male_population, individual_index)

                case 'Infection / Recovery': 
                    male_population = infection_function(male_population, individual_index)
        reactions.append(f'{sex}, {reaction_type}')


        t += tau
        iteration += 1
        
        time_series.append(t)
        herm_sizes.append(len(herm_population))
        male_sizes.append(len(male_population))

        alpha = parasite_genotype_frequencies(
            herm_population,
            male_population,
            kwargs['genotype_dictionary']
        )

        s_frequency , c_frequency = mating_type_frequency(
            s_frequency, 
            c_frequency, 
            herm_population, 
            male_population
        )
            
    return (herm_population,
            male_population,
            time_series,
            herm_sizes,
            male_sizes,
            s_frequency,
            c_frequency,
            reactions)
    

#Population Growth Parameters
REACTION_MAP = {
    "herm" : {
            0: "Death",
            1: "Birth",
            2: "Infection / Recovery"
        },
    "male" : {
        0: "Death",
        1: "Infection / Recovery"
    }
}

au = 0.006      #density dependance of uninfected individuals
ai = 1 * au   #density dependance of infected individuals

b0i = 1         #maximum birth rate of infected
b0u = 1 * b0i   #maximum birth rate of uninfected
d = 0.3         #death rate

ro = 0.5        #average fraction of males in a herm's offspring
rs = 0.9        #average fraction of males in a female's offspring

r = 0.3        #recombination rate

o_r = 0.5 

infection_effect_on_male = 0

T = 300
IP = 50

beta = 0.6
alpha = [random.randrange(0, 500) for _ in range(45)]
ssum = sum(alpha)
alpha = [alpha[z] / ssum for z in range (45)]


mp = population_dataframe_generator(40)
hp = population_dataframe_generator(40)
final_s_frequency = []
genotype_dictionary = two_locus_dictionary_function()

final_s_frequency = []
scheduled_events = [500]

herm_pop, male_pop, times, herm_sizes, male_size, s_freq, c_freq, reactions = \
    simulate_population(
        hp = hp, 
        mp = mp,
        genotype_dictionary = genotype_dictionary,
        t_max = 60, 
        scheduled_events = scheduled_events, 
        b0i = b0i,
        b0u = b0u, 
        ai = ai, 
        au = au,
        o_r = o_r, 
        rs = rs,
        ro = ro, 
        r = r, 
        beta = beta, 
        infection_effect_on_male = infection_effect_on_male, 
        d = d,
        alpha = alpha
        )
