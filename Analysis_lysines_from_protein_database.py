####################################################################################################################################################################
# ANALYSIS OF THE AMINO ACIDS COMPOSITION OF THE ENVIRONMENT OF LYSINES LISTED IN A PROTEIN DATABASE UNDER FASTA FORMAT(Version 1.0)
# CREATORS: LAB BOISVERT, UNIVERSITY OF SHERBROOKE (QC, CANADA)
# This work is licensed under Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) - 2024
####################################################################################################################################################################


# FORMATING THE DATABASE FILE (FASTA FORMAT) INTO A TEXT FILE FOLLOWING THE RULE "1 LINE = 1 PROTEIN"
#----------------------------------------------------------------------------------------------------------------------------------------------------------
print("Database cleaning in progress...")                   # Visual clue for the user to see which part of the program is running.
File = open("Uniprot_database_FASTA.txt","r")                             # Open the file containing the sequences of interest (it should be place in the same directory than the program file).
Raw_data = File.readlines()                                
File.close()                                                

Headers=[">"]                                               # Identify FASTA headers.
Numbers_of_proteins=0                                       # Count the number of proteins in the file.
Clean_File = open("Cleaned_database.txt","w")               
for line in Raw_data:                           
    if not any(Header in line for Header in Headers):       
        line=line.strip()                                   # Transform the file to obtain 1 line = 1 protein.                      
        Clean_File.write(line) 
    else:
        Clean_File.write("\n")                              # If a FASTA header is detected, a new line starts. 
        Numbers_of_proteins=Numbers_of_proteins+1           # 1 new line corresponds to a new protein.
print("Cleaning of the database is a success ! Number of proteins detected : ", Numbers_of_proteins)

Clean_File.close()                                          



#ISOLATION OF EACH LYSINE'S ENVIRONMENT (WINDOW GOING FROM -15 TO + 15 RESIDUES) 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
print("Isolating lysines ...")                              # Visual clue for the user to see which part of the program is running.
Initial_File=open("Cleaned_database.txt","r")               # Re-open the file previously cleaned but in a "reading" mode to protect those new data.
Initial_Database = Initial_File.readlines()                 
Modification_windows=open("Lysine_environment_raw.txt","w")    

Numbers_of_windows=0                                        # Count the number of windows (equal to the number of lysine) in the cleaned database.

for Database_sqce in Initial_Database:                      
    l=len(Database_sqce)                                    # Verify the length of each line (ie each sequence).
    for i in range (l):                                     
        if Database_sqce[i]=="K":
            if i<15:                                        # If the lysine is to close to the N-term end.  
                Modification_windows.write((15-i)*"X"+Database_sqce[0:i+16]+"\n")       # Add the letter "X" for missing residue.        
            elif l-i<=16:                                   # If the lysine is to close to the C-term end.
                e=15-(l-i)+2                                # Determine the number of "X" to add in order to reach a sequence of 31 residues.
                Modification_windows.write(Database_sqce[i-15:l-1]+ e*"X"+"\n")
            else:                                           
                Modification_windows.write(Database_sqce[i-15:i+16]+"\n")       
            Numbers_of_windows=Numbers_of_windows+1             
        
#print(Numbers_of_windows)

Initial_File.close()
Modification_windows.close()


Modification_windows_2=open("Lysine_environment_raw.txt","r")

Modification_windows_FINAL=open("Isolated_sequences.txt","w")    
Rubbish=0                                                   

Modification_windows_read=Modification_windows_2.readlines()

for Windows_line in Modification_windows_read:                   # Get rid of empty lines that may have appeared during formating the original file
    if len(Windows_line)<=30:

        Rubbish=Rubbish+1
    else:
        Modification_windows_FINAL.write(Windows_line)

#print("Number of empty lines : ", Rubbish)
print("Each lysine's environment was successfully isolated. Starting environment analysis.")


Modification_windows_FINAL.close()



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
File_Sequences = open("Isolated_sequences.txt","r")        # Open the document containing the sequences of interest (Cleaned sequences contening only 31 residues).
Raw_data2 = File_Sequences.readlines()


Clean_raw_data2 = []                     # New matrix for used to clean isolated sequences.              
for line in Raw_data2:                   # Remove "\n" from the original sequence.
    Clean_raw_data2.append(line.strip())

Sqce=[list(individual) for individual in Clean_raw_data2]    
Total_sites=len(Clean_raw_data2)    
print("Analyzing "+str(Total_sites)+" sites...") # Visual clue for the user which indicates the number of sites detected by the program.

Apparition_matrix = [[0 for i in range (21)] for j in range (31)]          # Initializing the matrix counting the number of apparition for each position and each residue (20 amino acids + gaps).

for y in range(Total_sites):             
    for z in range (30):                
        j=check_aa(Sqce[y][z])          # Identification of the residue.
        Apparition_matrix[z][j]=Apparition_matrix[z][j]+1                # Incrementation.
            
print(Apparition_matrix)                       # Return the final apparition matrix.

File_Sequences.close()
