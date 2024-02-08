####################################################################################################################################################################
# ANALYSIS OF THE AMINO ACIDS COMPOSITION OF RESIDUES SURROUNDING MODIFIED LYSINES IDENTIFIED IN DIGLY IMMUNOPRECIPITATION (Version 1.0)
# CREATORS: LAB BOISVERT, UNIVERSITY OF SHERBROOKE (QC, CANADA)
# This work is licensed under Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) - 2024
####################################################################################################################################################################


#DEFINE A FUNCTION WHICH IDENTIFY THE AMINO ACID X IN THE CHARACTER CHAIN
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
def check_aa(x):                                         
    if x=="A":
        AA=0
    elif x=="C":
        AA=1
    elif x=="D":
        AA=2
    elif x=="E":
        AA=3
    elif x=="F":
        AA=4
    elif x=="G":
        AA=5
    elif x=="H":
        AA=6
    elif x=="I":
        AA=7
    elif x=="K":
        AA=8
    elif x=="L":
        AA=9
    elif x=="M":
        AA=10
    elif x=="N":
        AA=11
    elif x=="P":
        AA=12
    elif x=="Q":
        AA=13
    elif x=="R":
        AA=14
    elif x=="S":
        AA=15
    elif x=="T":
        AA=16
    elif x=="V":
        AA=17
    elif x=="W":
        AA=18
    elif x=="Y":
        AA=19
    else:                               # If no amino acid is read, return an additional value indicating a gap.
        AA=20
    return(AA)                          # Return the value corresponding to the amino acid.


                                        
# GENERATION OF THE MATRIX INDICATING THE NUMBER OD AMINO ACIDS IDENTIFIED FOR EACH POSITION
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
File_Sequences = open("Isolated_sequences_diGly_IP.txt","r")        # Open the document containing the sequences of interest (Cleaned sequences contening only 31 residues).
Raw_data2 = File_Sequences.readlines()


Clean_raw_data2 = []                     # New matrix for used to clean isolated sequences.              
for line in Raw_data2:                   # Remove "\n" from the original sequence.
    Clean_raw_data2.append(line.strip())

Sqce=[list(individual) for individual in Clean_raw_data2]    
Total_sites=len(Clean_raw_data2)    
print("Analyzing "+str(Total_sites)+" sites...") # Visual clue for the user which indicates the number of sites detected by the program.

Apparition_matrix = [[0 for i in range (21)] for j in range (31)]          # Initializing the matrix counting the number of apparition for each position and each residue (20 amino acids + gaps).

for y in range(Total_sites):             
    for z in range (31):                
        j=check_aa(Sqce[y][z])          # Identification of the residue.
        Apparition_matrix[z][j]=Apparition_matrix[z][j]+1                # Incrementation.
            
print(Apparition_matrix)                       # Return the final apparition matrix.

File_Sequences.close()
