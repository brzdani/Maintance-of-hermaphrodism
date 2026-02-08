import numpy as np
import pandas as pd 
import random
import math 
from itertools import product
from itertools import combinations_with_replacement
import matplotlib.pyplot as plt
from tqdm import tqdm


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


def  two_locus_index_dictionary_function(allele_number = 3):
    counter = 0
    genotypes = []
    values = []
    for i in range(allele_number):
        for j in range(allele_number):
            for z in range(allele_number):
                for y in range(allele_number):
                    if i == z and j == y:
                        genotypes.append(f"A{i+1} B{j+1} / A{z+1} B{y+1}")
                        values.append(counter)
                        counter += 1
                    else:
                        if genotypes.count(f"A{i+1} B{j+1} / A{z+1} B{y+1}") == 0:
                            genotypes.append(f"A{i+1} B{j+1} / A{z+1} B{y+1}")
                            values.append(counter)
                            counter += 1
    index_dictionary = {values[i]: f'{genotypes[i]}' for i in range(len(genotypes))}
    return index_dictionary

def genotype_converter (genotype, genotype_dictionary):
    '''
    Parametes 
    ---------
    Genotype: a string containing the individual's genotype

    
    Return
    ------
    Returns the index of each genotype of two loci with 3 alleles. (Ranging from 0 to 44)  
    '''

    return genotype_dictionary.get(genotype, 1000)

def index_to_genotype (x, index_dictionary):
    '''
    Parametes
    ---------
    ّIndex: a string containing the genotype's index (Ranging from 0 to 44) 

    
    Return
    ------
    Returns the genotype related to the given index.  
    '''

    return index_dictionary.get(x, 10000)

def parasite_genotype_frequencies(herm_population, male_population, iteration = None):


    population = pd.concat((herm_population, male_population), ignore_index=True)
    infected = population[population['Parasite status'] == 'infected']

    if len(infected) == 0: 
        genotype_dic = two_locus_dictionary_function()
        return [0.0] * len(set(genotype_dic.values()))
    
    genotype_dic = two_locus_dictionary_function()
    alpha = [0.0] * len(set(genotype_dic.values()))

    for genotype in infected['Genotype']: 
        u = genotype_converter(genotype, genotype_dic)
        alpha[u] += 1/len(infected)

    return alpha

def infection_function (population, individual_index):
    '''
    Parameters
    ----------
    A data frame called population contianing all details about the population 
    (herm or male) on which we want to apply the infection proccess. Also gets 
    the index regarding the individual that their infection status will be changed.  


    Return
    ------
    New data frame containing the new population informations. 
    '''

    population.at[individual_index , 'Parasite status'] = 'uninfected' \
        if population.at[individual_index, 'Parasite status'] == 'infected'\
            else 'infected'
    return population

def death_function (population, individual_index):

    '''
    Parameters
    ----------
    A data frame called population contianing all details about the population (herm or male)
    on which we want to apply the death proccess. Also gets the index regarding the individual 
    that will be removed from the population via death function.  


    Return
    ------
    New data frame containing the new population informations. 
    '''

    population = population.drop(individual_index).reset_index(drop = True)
    population['Index'] = population.index
    return population

def mating_male_index (male_population, infection_effect_on_male): 
    ''' 
    Parameters
    ----------
    gets the male population. 


    Return
    ------
    Return which male will mate with the hermaphrodite based on the male infection status and parasite effect on males. 
    '''
    weights = []
    for status in male_population['Parasite status']:
        if status == 'infected': 
            weights.append(1-infection_effect_on_male)
        else:
            weights.append(1.0)
    weights = np.array(weights)
    weights /= weights.sum()

    return np.random.choice(len(male_population), p=weights)



def recombination (gen, r):
    
    ''' 
    Decides if recombination happens or not. 
    Return a list of alleles on the chromoses (indeces 0, 1 on the first chromosome
    and indeces 2 , 3 on the second chromosome.)
    '''
    if isinstance(gen, str):
        gen = gen.split()

    if np.random.rand() < r:
        gen[1], gen[4] = gen[4], gen[1]

    return gen

def newborn_gender_generator (r):
    '''        
    Decides what will bethe gender of the neworn based on sex ratio. 
    ''' 

    return 'herm' if np.random.rand() < r else 'male'

def newborn_genotype_generator(herm_gen, male_gen, r):
    '''   
    Parameters 
    ----------
    Gets genotype of male and female (producers of sperm and ovule are called respectively male and female). 

    
    Return
    ------
    Newborn genotype. 
    '''
    herm_gen = recombination(herm_gen, r)
    male_gen = recombination(male_gen, r)

    sperm = np.random.randint(2)
    ovule = np.random.randint(2)

    sperm_genotype = " ".join(male_gen[2 * sperm : 2 * sperm + 2])
    ovule_genotype = " ".join(herm_gen[2 * ovule : 2 * ovule + 2])

    return f'{sperm_genotype} / {ovule_genotype}'

def mating_type_generator (herm_population,
                           individual_index,
                           male_population,
                           index_of_mating_male):
    
    mt_ovule = herm_population.at[individual_index, 'Mating type'].split()
    mt_ovule = mt_ovule[2*np.random.randint(2)]

    mt_sperm = male_population.at[index_of_mating_male, 'Mating type'].split()
    mt_sperm = mt_sperm[2*np.random.randint(2)]

    

    return f'mt_sperm / mt_ovule'


def birth_function (herm_population, male_population, individual_index, r,
                    infection_effect_on_male = 0, ro = 0.5, rs = 0.9, o_r = 0.5):

    male_index = 0
    newborn_parasite_status = 'uninfected'
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
        male_gen = herm_population.at[individual_index, 'Genotype']
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
            len(male_population), 
            newborn_genotype, 
            mt, 
            newborn_parasite_status
        ]

    if newborn_gender == 'herm':
        herm_population.loc[len(herm_population)] = [
            len(herm_population),
            newborn_genotype,
            mt,
            newborn_parasite_status
        ]

    return herm_population, male_population
    
def genotype_converter (genotype, genotype_dictionary):
    '''
    Parametes 
    ---------
    Genotype: a string containing the individual's genotype

    
    Return
    ------
    Returns the index of each genotype of two loci with 3 alleles. (Ranging from 0 to 44)  
    '''

    return genotype_dictionary.get(genotype, 1000)

def mating_type_frequency(s_frequency, c_frequency, herm_population, male_population, iteration):

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

    s_frequency.append(s_count / N)
    c_frequency.append(c_count / N)

    return s_frequency, c_frequency  


def genotype_generator (): 
    '''
    Parameters
    ----------
    None. 

    Returns: 
    -------
    Generates a random genotype between all possible genotypes in case of 2 loci with 3 alleles 
    '''
    genotype_dic = two_locus_index_dictionary_function()
    index = random.randint(0, len(genotype_dic))
    genotype = index_to_genotype(index, genotype_dic)

    return genotype

def population_dataframe_generator(IP): 
    return pd.DataFrame({
        "Genotype": [genotype_generator() for _ in range(IP)],
        "Mating Type": ["C / C"] * IP, 
        "Parasite status": np.random.choice(
            ["infected", "uninfected"],
            size = IP
            )
    })

def tau_generator(a): 
    '''
    create a random time till the next event based on the reaction rate a 
    '''
    if a <= 0:
        return math.inf
    return np.random.exponential(1/a)


def reaction_type_converter(sex, i):
    return REACTION_MAP[sex].get(i, "Unkown") 



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



for n in tqdm(range(0, 30)):
    herm_population = hp.copy(deep=True)
    male_population = mp.copy(deep=True)

    male = [len(male_population)]  
    herm = [len(herm_population)] 

    s_frequency = list()
    c_frequency = list()
    s_frequency , c_frequency = mating_type_frequency(s_frequency, c_frequency, hp, mp, 0)

    t = 0
    iteration = 0
    time = [t]
    reactions = [0]
    while t < 60:
        if iteration == 500: 
            herm_population.at[0, 'Mating type'] = "S / S"

        ta_herm = list()
        ta_male = list()
        for i in herm_population['Index']:

            
            if herm_population.at[i,'Parasite status'] == 'infected':
                ta_herm.append(tau_generator(d*2))
                ta_herm.append(tau_generator(b0i - ai * (herm[iteration] + male[iteration])))
                ta_herm.append(tau_generator(beta))
                
            if herm_population.at[i,'Parasite status'] == 'uninfected':
                ta_herm.append(tau_generator(d))
                ta_herm.append(tau_generator(b0u - au * (herm[iteration] + male[iteration])))
                ta_herm.append(tau_generator(alpha[genotype_converter(herm_population.at[i,'Genotype'], genotype_dictionary)]))
                
        for z in male_population['Index']:
            
            if male_population.at[z,'Parasite status'] == 'uninfected':
                ta_male.append(tau_generator(d))
                ta_male.append(tau_generator(alpha[genotype_converter(male_population.at[z,'Genotype'], genotype_dictionary)]))
            if male_population.at[z,'Parasite status'] == 'infected':
                ta_male.append(tau_generator(d*2))
                ta_male.append(tau_generator(beta))

        sex_reaction = 'Herm' if min(ta_herm) <= min(ta_male) else 'Male'
        
        if sex_reaction == 'Herm':
            tau = min(ta_herm)
            reaction_type = reaction_type_converter('herm', ta_herm.index(min(ta_herm)) % 3)
            individual_index = ta_herm.index(min(ta_herm)) // 3
            if (reaction_type == 'Birth'):
                reactions.append(f'{sex_reaction} , {reaction_type}')
                herm_population, male_population = birth_function(herm_population, male_population,
                                                                   individual_index, r, infection_effect_on_male, ro, rs, o_r)
            if (reaction_type == 'Death'):
                reactions.append(f'{sex_reaction} , {reaction_type}')
                herm_population = death_function(herm_population, individual_index)
            if (reaction_type == 'Infection / Recovery'):
                reactions.append(f'{sex_reaction} , {reaction_type}')
                herm_population = infection_function(herm_population, individual_index)
        elif sex_reaction == 'Male':
            tau = min(ta_male)
            reaction_type = reaction_type_converter('male', ta_male.index(min(ta_male)) % 2)
            individual_index = ta_male.index(min(ta_male)) // 2
            if (reaction_type == 'Death'):
                reactions.append(f'{sex_reaction} , {reaction_type}')
                male_population = death_function(male_population, individual_index)
            if (reaction_type == 'Infection / Recovery'):
                reactions.append(f'{sex_reaction} , {reaction_type}')
                male_population = infection_function(male_population, individual_index)

        male.append(len(male_population))
        herm.append(len(herm_population))
        

        t = t + tau
        time.append(t)
        iteration += 1

        alpha = parasite_genotype_frequencies(herm_population, male_population, iteration) 
        s_frequency , c_frequency = mating_type_frequency(s_frequency, c_frequency, herm_population, male_population, iteration)
    
    final_s_frequency.append(s_frequency[-1])

    xpoints = np.array(time)
    ycpoints = np.array(c_frequency)
    yspoints = np.array(s_frequency)
    ysplot = plt.plot(xpoints, yspoints, color='m', alpha= 0.7, linewidth = 0.7)



plt.title('Population Size')
plt.xlabel('Time')
plt.ylabel('Relative Frequency')
plt.grid(alpha = 0.3)
plt.show()