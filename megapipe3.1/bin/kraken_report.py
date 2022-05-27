import os

# looking through directory and getting names of all folers and files
# directly in that directory

# path to directory
input_path = "/n/data1/hms/dbmi/farhat/rollingDB/cryptic_output/"

# makes a list
files = os.listdir(input_path)
ids = []
no_kraken = []
for id in files:

    if id[0]== 'E' and (os.path.isdir(input_path + str(id) + '/kraken')):
        ids.append(id)
    else:
        no_kraken.append(id)

# for counting each row 
counter = 0
# for counting each good sequence that gets classified to 1773
sequence = 0

# counting number of strains that have 90% and above good sequences
good_runs = 0

# tracking strains that didnt pass
bad_ids = []

# opening each kraken file
for id in ids:
    sequence = 0
    counter = 1
    with open(input_path + str(id) + '/kraken/'+ str(id) +'.krkn', 'r') as file:
        for line in file:
            row = line.split('\t')
            counter += 1
           
            # lazy way of making sure the row at least continas the third column
            if len(row) >= 3:
                if row[2] == '1773':
                     sequence += 1
        # checking if more than 90% of sequences are good             
        if (sequence / counter) >= 0.9:
             good_runs += 1
        # making a list of the ids that didnt pass
        else:
            bad_ids.append(id)
   
print(f'Number of succesfully classified isolates: {good_runs}')
print(f'Number of failers: {len(bad_ids)}')


with open('cryptic_kraken_report.txt', 'w') as report:
    report.write('Good runs: ' + str(good_runs))
    report.write('\n')
    report.write('Isolates that didnt classifiy 90% and above: ')
    for id in bad_ids:
        report.write(id)




