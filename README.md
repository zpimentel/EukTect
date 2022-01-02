# EukTect
Development of a SVM to identify eukaryotic contigs from bacterial, archaeal, and viral contigs. This implementation was run on a laptop and thus does not deeply sample many genomes - likely limiting its accuracy.

## Step 1. Prepare the environment
```
$ git clone git@github.com:zpimentel/EukTect.git
$ conda env create -f environment.yml
$ source activate EukTect
```

## Step 2. Train and test the model
This step randomly selects and downloads 50 eukaryotic genomes, 20 viral genomes, 15 bacterial genomes, and 15 archaeal genomes to use for training of a SVM. I go through up to 50 kb of each genome, breaking them into 5 kb long non-overlapping chunks/contigs. Then, I use a sliding window to count canonical kmers (word size of 5) in each contig - assigning each contig as eukaryotic or other (meaning bacterial, archaeal, or viral) in the process. Then, I calculate the frequencies of the canonical kmers and use this to train a SVM. This step uses the NCBI RefSeq summary table (https://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt) downloaded on 2022-01-01. This table is included in the git repo because it changes periodically and it is randomly sampled; thus, different results could be obtained if different versions of the table were used.
```
$ jupyter lab
$ python 01_train_and_test.py
```
I reached a fairly high accuracy in the testing set (96.7%). However, it must be remembered that while the testing set contigs were not seen by the model, they come from the same genomes that the training contigs did. Therefore, there are very likely to be biases associated with this testing data set.  

## Step 3. Validation on genomes not seen by the model
This steps performs validation of the model using new genomes that the model has not seen before. I used 10 bacterial and 10 eukaryotic genomes. While this is not many, I created thousands of contigs that were independently tested from each of these genomes. Therefore, thousands of contigs were ultimately tested. Another approach would be to download more genomes and to sample a smaller number of contigs from each one.
```
$ jupyter lab
# Open and run 02_validation.ipynb
```
From this test, I reached an accuracy of 100%. However, I know there is still work to be done such as inclusion of some GC-rich bacteria in the training data.   

## Step 4. Test the model on a user-defined sequence
This step allows a user to input a sequence and a prediction will be generated whether it is of eukaryotic origin. Sequence should generally be > ~ 3 kb. Currently unknown how archaeal or viral contigs would be classified.
```
$ jupyter lab
# Open and run 03_apply_model.ipynb
```

