"""

Author: Hui-Heng Lin, Phd. Https://orcid.org/0000-0003-4060-7336




Purpose (Network link prediction via machine learning): To create the drug-target interaction pairs' machine learning matrix data (feature vectorsets), and to use machine learning method (Random forest) to predict the probable drug-target interaction links. 



Notes: 
(1) Before starting, one needs to prepare a list drug-target interaction pairs by yourself.  

(2) For negative reference instances, one can consider a drug and one of its non-interactive protein as an negative training/validating instance. Through random sampling, multiple such kind of negative (or non-interactive drug-target ) pairs could be obtained.

(3) As matrice would be combined, it is of vital importance to match the correct order of drug and its corresponding interative target protein in the item list, i.e., the .SDF file for drug data and the .FASTA file for protein data,  respectively. Otherwise, the orders and pairs will be wrong and it certainly leads to wrong results. 

"""




### installation of required packs ####
install.packages("BiocManager")
BiocManager::install("ChemmineR") # required for analyzing drug data
BiocManager::install("protr") # required for analyzing drug target(protein) data


# load the drug dataset first, and transfer them into 1024 vectors
library(ChemmineR) # load the package required for processing drug dataset


drugs<-read.SDFset("N:/your_local_path/drug_structure_data_file.sdf") # need to collect .sdf (drug molecular structure data) file so as to load the drug structure data.  Note that multiple molecules' structure data could be saved in a single .SDF file.


drugs # examine if the drug molecules were loaded properly
view(drugs[1]) # same above


drugfp<-desc2fp(sdf2ap(drugs)) # Initially, convert the drug SDF data into atompair descriptor data via sdf2ap( ); Next, use desc2fp( ) function to convert atompair descriptor data into fingerprint data.


drugfp@fpma # check and access the matrix of fingerprint data
class(drugfp@fpma) # it was a matrix
dim(drugfp@fpma) # checking if the data were correct via showing the dimension of the fingerprint matrix data of drug molecules



# Attach a columns of label data to above matrix for machine learning purposes. E.g., the 0/1 so as to indicate the yes/no or other types of binary classification. It could be done via the column-binding function --- cbind( ) 


# create a 1-colunm matrix data first
lab<-matrix(data= NA, X,1) # X should be the number of your instance for machine learning 
dim(lab) ; head(lab); tail(lab) # confirm the dimension of the matrix, and the elements inside 



m2<-cbind(lab,drugfp@fpma) # notice that both variables here were matrix type data.
dim(cbind(lab,drugfp@fpma)); head(m2); tail(m2) # check the correctness of the dimension of matrix and also part of the data



# Analyzing protein dataset

library(protr) # load package required for analyzing the protein dataset

# load the protein sequence dataset first
tg<-readFASTA("F:/your_local_path/protein_seq.fasta") # protein sequences are usually stored/saved as .FASTA format file.

class(tg) # it was a list type variable!

length(tg); head(tg); tail(tg) # checking the connectness of number of protein sequences, and parts of protein data loaded. Notice that,  headers of .FASTA format: ">abcd_xxx" will be lost! Therefore, as mentioned at the beginning notes, it is of vital importance to mark or adjust the correctness of the order of the protein sequences in advance. Otherwise, it will be difficult to identify whether your drug-target interaction pairs are correct or not.

tg[1]; tg[2];  tg[[1]]; tg[[2]] # try to access the elements and check
 # hehe the same outputed as [1] o!



# Don't forget to use below function to check whether the protein seq are in the 20 default types of amino acids. If unrecognized alphabetical letter (amino acid other than the 20 typical types), errors will happen during computations. 
protcheck(x) # x is a character vector, as the input protein sequence


# to obtain numeric representations of protein seqeuences
extractAAC(tg[[1]]) #  This is the function for getting the 20-type amino acid composition of a protein sequence. I used it as the only example function in this tutorial.

""" More advanced functions include but not limited to, e.g., extractGeary( ) , extractAPAAC( ), extractCTDC( ), etc. For more functions and details, please refer to the manual of protr package.

One can use multiple numeric representation functions for protein sequences so as to generate different vectors of a protein sequence. And then can use cbind( ) function to integrate multiple vectorsets together for machine learning purpose. Hereby one can customize your featureset of machine learning, and better model might be obtained via optimization of featureset. 

"""



# we could use loop to vectorize multiple protein sequence efficiently.

protmtr<-matrix(NA,x,20); # an empty matrix was initially created for storing vectorset in advance. The "x" was the number your instance for machine learning 

# the loop
for (i in 1:x){ protmtr[i,]<- extractAAC(tg[[i]]);  print(paste(i,"--th round, ends----"))}  # The "x" was the number your instance for machine learning 



# check if the resultant matrix data were sound and fine
head(protmtr) 



# to combine all the matrix into full  matrix for machine learning !
fulmat<-cbind(m2,protmtr) ; dim(cbind(m2,protmtr)) ; dim(fulmat) # combination and checking


# install and use random forest 
install.packages("randomForest") # installation 
library(randomForest) # load the pack


# train a model
rf<-randomForest(x = fulmat[,x:y], # Submit all the instances to train a model. x and y were two column numbers indicating the range of feature selected for training model. 
                 y = as.factor(fulmat[,z]) ) # z is the number of the labelset column to be specified.  Note that if as.factor( ) did not present here, it would become regression analysis instead of machine learning analysis.


# Above was the training of model without manually specified validation set input, in which the function will performance cross-validation itself. Nevertheless, in fact, the randomForest( ) function allows users to manully specify or provide their desired validation set. See example below:
rf_v<-randomForest(x = fulmat[x:y,z:a], # x and y here were the number/range of instances including only the training set while not covering the validation set . And  "z:a"  was the range of specified featureset for training 
                 y = as.factor(fulmat[x:y,k]), # k was the column number of the machine learning labelset for training  and model generation 
                 xtest=fulmat[u:v,z:a], # u and v were the number or range of instances for model performance validating purpose, u:v did not include the range of trainingset
                 ytest=as.factor(fulmat[u:v,k]), ntree=500 )  # Here only list basic arguments of randomForest( ). Further and more advanced arguments could be found at the manual of randomForest pack


"""
The running of above code will return information like:






Call:
 randomForest(x = fulmat[, 2:1045], y = as.factor(fulmat[, 1]),      ntree = 500) 
               Type of random forest: classification
                     Number of trees: 500
No. of variables tried at each split: 32

        OOB estimate of  error rate: 14.29%

Confusion matrix:
   0  1 class.error
0 81 10   0.1098901
1 16 75   0.1758242





Above confusion matrix helps us to evaluate the performance of the model. Amongst, value 81 was in the position of 00 (number of true negative instance), indicating there were 81 true negative instances identified by the model. Similarly, 10 was the number of the false positive instances, 16 was number of false negative instances, and 75 was the number of true positive recognized by the model. 

So the discovery rates of true positive, true negative, false positive and false negative, could be further calculated from these values.

For example, the true positive rate = number of true positive instances / total number of positive instances = 75 / (75 + number of false negative instances) = 75 / (75+16) = 0.82417

For further instructions and programming codes about how to convert confusion matrix into discovery rates and teh value of area under receiver-operating characteristic curves (AUROC), see here  https://github.com/devgithwaar/compute_true_positive_rates_AUCs   


"""



# PS : Variables in R environment could be exported and saved in local computer file. Sometimes this helps a lot and can greatly reduce the effort for various purposes like reloading data or restoring models. Belows are several examples for exporting data into local computer files.

# When .FASTA file was modified in R, you may want to export and save the changed one into computer file for next time's convenient usage.

# file Connection way
fc<-file("N:/local_computer_path/output.fasta"); writeLines(as.character(tg),fc) ; close(fc) #  the protein sequences were output and saparated by lines! BUT notice that headers "> xxxx" of sequence info were lost through such operation




# matrix varaibles could be exported as the .csv (comma separated value) files in local computer pathway. 

# for instance, to  export the drug 1024 fingerprint matrix combined with label 1/0 column data --- can use write.csv( ) directly!
write.csv(m2, "N:/local_path/drugFPlabelMatrix.csv", sep=",") 

# to export all the protein seqs' amino acids composition matrix 
write.csv(protmtr, "N:/local_path/proteinProtRAACMatrix.csv") # 


# also export the full matrix for machine learning to local file 
write.csv(fulmat, "N:/local_path/fullMatrix.csv")
