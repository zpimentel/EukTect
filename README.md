# EukTect
Development of a SVM to classify bacterial and eukaryotic contigs using canonical kmers.

## Step 1. Prepare the environment
```
$ git clone git@github.com:zpimentel/EukTect.git
$ conda env create -f environment.yml
$ source activate EukTect
```

## Step 2. Train and test the model
This step randomly selects 50 bacterial and 50 eukaryotic genomes to use for training of a SVM. I go through up to 1 Mb of each genome, breaking them into 50 kb long non-overlapping chunks/contigs. Then, I use a sliding window to count canonical kmers (word size of 5) in each contig - assigning each contig as eukaryotic or bacterial in the process. Then, I calculate the frequencies of the kmers and use this to train a SVM.
```
$ jupyter lab
# Open and run 01_train_and_test.ipynb
```
I reached a very high accuracy in the testing set (99.6%). However, it must be remembered that while the testing set contigs were not seen by the model, they come from the same genomes that the training contigs did. Therefore, there are very likely to be biases associated with this testing data set.  

Things to improve: Rather than using randomly selected eukaryotic and bacterial genomes for training select a set of taxonomically diverse genomes, make more efforts to make sure there is an even amount of bacterial and euk data in the training data. Specifically, I know the model performs poorly on GC-rich bacterial genomes which weren't seen in the training set, so inclusion of some of those could help.

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

