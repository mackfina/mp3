import os

# looking through directory and getting names of all folers and files 
# directly in that directory

# path to directory
input_path = "/n/data1/hms/dbmi/farhat/fastq_db/MIC-ML_Consortium/Wadsworth"

# makes a list
files = os.listdir(input_path)

test_list = ['ERR2510808', 'ERR2199868']
# make a file with those ids

# location for the id's.txt file
output_path = "/n/data1/hms/dbmi/farhat/rollingDB/MIC-ML_output/Wadsworth/logs/run1/ids.txt"
#counter for only getting i amount of strains for batching
i = 0

# removing extra stuff from end:

ids = []

for id in files:
    ids.append(id[:-16])
# removing duplicate ids:
unique_ids = []
for id in ids:
    if id not in unique_ids:
        unique_ids.append(id)



# write ids from directory to output
with open(output_path, 'w') as f:   
        for id in unique_ids:
            f.write(id+ '\n')
        
for item in unique_ids:
    print(item)
print(len(unique_ids))
print("Done gettingIDs!")


