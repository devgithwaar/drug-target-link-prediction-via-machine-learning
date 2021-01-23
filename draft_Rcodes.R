"""
7th december 2020@shaoguan yuebei hospital 3-13f-1341 raum 
Rstudio win10 hwaarwhy  [updated on 8th december 2020@shaguan shiwei, win8.1 hp omen 17's Rstudio ! ]

Author Huiheng lin phd

Purpose: create the drug-target interaction pairs' machine learning matrix data(feature vectorsets)


"""


## 
### ## !!!! odel DO I still need to export the random forest model(s)? into local file ????????? shit, it was in the HP omen 17, not in other PC!! need to retrain and load many data if want to get the model here !
###  or I should save R.data and make it portable ?? store in USB disk or mobile harddisk?





### installation of required packs ####
install.packages("BiocManager")
# and then
BiocManager::install("ChemmineR")
BiocManager::install("protr") # read document to see what descriptor functions are available


# load the drug set first, and transfer them into 1024 vectors

drugs<-read.SDFset("N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/drugSDFdataset/IntegratedAlldrugsAntiviralSDFsSortedAsDrug-TargetPairs.sdf")

drugs # cool 182 molecules were loaded

view(drugs[1]) # cool. checked

drugfp<-desc2fp(sdf2ap(drugs)) # cool. 

drugfp@fpma
dim(drugfp@fpma) # cool 182 x 1024 matrix!

head(drugfp@fpma[,])
head(drugfp@fpma[1,])
drugfp@fpma[1:8,] # odel show numbers, not 0/1
dim(drugfp@fpma[1:8,]) #  8 x 1024, seemed normal

drugfp@fpma[1:8,3] #  0/1 are shown
col(drugfp) # error!


# attach a label columns for machine learning, i.e, the 0/1 or yes/no. using  column bind, cbind( ) function.

""" machine learning data format of the excel file:
[1:73] are the positive training drug-target pairs --- should be labelled as 1

[74:146] are the negative training drug-target pairs -- should be labelled as 0!

[147:164] are the positive validating set --- should be labelled as 1

while [165:182 ] are the negative validating set of drug-target pairs(unlinked pairs randomly sampled , hehe) ---- should be labelled as 0 !

"""

# create the 1-colun first
lab<-matrix(data= NA, 182,1)
dim(lab) # confirmed, 182 rows x 1 column
head(lab) # hehe all NA
tail(lab) # NAs, too

lab[1:73,]<-1; lab[1:73,] # cool all be came 1!
lab[74:146,]<-0; lab[73:146,] # cool all be came 0!
lab[147:164,]<-1; lab[147:164,] # cool
lab[165:182,]<-0; lab[165:182,] # cool

cbind(lab,drugfp) # odel. bound, but show ???? question marks!!
cbind(lab,drugfp@fpma) # odel. still not showing 0/1 just large umbers!
dim(cbind(lab,drugfp@fpma)) # 182 x 1205 , seemingly fine
cbind(lab,drugfp@fpma)[1:70,3:10] # cool! showed 0/1!
cbind(lab,drugfp@fpma)[70:100,1:5] # cool! showed 0/1!



library(protr)
tg<-readFASTA("F:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/182orderedProteinTargets4MachineLearning5dec2020.fasta") # cool. loaded no error---------this command was written in server at hospital (win 2016 server) -----now in HP omen 17 the disk character changed from F to N! hehe1 

tg<-readFASTA("N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/182orderedProteinTargets4MachineLearning5dec2020.fasta") # cool.

length(tg) # 182 ,cool

head(tg) # shit headers all lost! shit, no headers can be read by this function!! some words were read, shit cannot recognize full?

tail(tg)

tg[1]
tg[2]
tg[[1]] # hehe the same outputed as [1] o!
tg[[2]] # hehe the same outputed as [2] o! seq only?



extractAAC(tg) # error. odel a character vector?

extractAAC(tg[1]) # error!
extractAAC(tg[[1]]) # cool! [[1]] seq only worked!

extractAAC(character(tg)) # error. invalid length?
extractAAC(as.character(tg)) # error


extractGeary(tg[[1]]) # heiyo outputed!


protmtr<-matrix(NA,182,20); 
for (i in 1:182){ protmtr[i,]<- extractAAC(tg[[i]]);  print(paste(i,"--th round, ends----"))}  # error in 31 st amino acid seqeuence ! X unrecognized amino acids presented!

tg[31] # odel. showed , but no X found!? the problematic one is 30-th?

tg[30] # also not found X, test by the function directly

extractAAC(tg[[30]]) # no problem
extractAAC(tg[[31]]) # shit, error ! X presented !----------shit!! no X ! other than X's amino acids presents?
extractAAC(character(tg[[31]])) # shit, still error !

# 6 strange letters excluded from 20 types of amino acids. B, J, O ,U , X, Z -----------shit checked tg[[31]] without these letter Masaka error is due to the  shit emtpy spaces!!?



extractAIC(character(tg[[31]])) # errror !
extractAIC(character(tg[[31]])) # errror !

extractAPAAC(tg[[31]]) # hehe error again, due to the X unrecognized.

# manual removing the space!
tg[[31]]
class(tg[[31]]) # character type

tg[[31]]<-"MEDFVRQCFNPMIVELAEKTMKEYGEDLKIETNKFAAICTHLEVCFMYSDFHFINEQGESIIVELGDPNALLKHRFEIIEGRDRTMAWTVVNSICNTTGAEKPKFLPDLYDYKENRFIEIGVTRREVHIYYLEKANKIKSEKTHIHIFSFTGEEMATKADYTLDEESRARIKTRLFTIRQEMASRGLWDSFRQSERGEETIEERFEITGTMRKLADQSLPPNFSSLENFRAYVDGFEPNGYIEGKLSQMSKEVNARIEPFLKTTPRPLRLPNGPPCSQRSKFLLMDALKLSIEDPSHEGEGIPLYDAIKCMRTFFGWKEPNVVKPHEKGINPNYLLSWKQVLAELQDIENEEKIPKTKNMKKTSQLKWALGENMAPEKVDFDDCKDVGDLKQYDSDEPELRSLASWIQNEFNKACELTDSSWIELDEIGEDVAPIEHIASMRRNYFTSEVSHCRATEYIMKGVYINTALLNASCAAMDDFQLIPMISKCRTKEGRRKTNLYGFIIKGRSHLRNDTDVVNFVSMEFSLTDPRLEPHKWEKYCVLEIGDMLIRSAIGQVSRPMFLYVRTNGTSKIKMKWGMEMRRCLLQSLQQIESMIEAESSVKEKDMTKEFFENKSETWPIGESPKGVEESSIGKVCRTLLAKSVFNSLYASPQLEGFSAESRKLLLIVQALRDNLEPGTFDLGGLYEAIEECLINDPWVLLNASWFNSFLTHALS"; tg[[31]] # the space is gone

extractAAC(tg[[31]]) # cool !!! shit fuck you!! really due to the shit space!!



protmtr<-matrix(NA,182,20); for (i in 1:182){ protmtr[i,]<- extractAAC(tg[[i]]);  print(paste(i,"--th round, ends----"))} # hehe encountered error in 72th!
tg[[72]] # haha, checked fuck you another space presented 

tg[[72]]<-"MSWAKQRVPFLDDDDGEEENDVQDDVDSPVPTRPLVIDEDAEPAAGTSGGLEGGGGDDEDGEDGHALPDLDDDLLLQFEPMLPRVYDLLLPSLDARLNFVNAGQKYAAFLKYVHGDCATCSHGEILREKTQLLTAIVSKLMDINGILEGKDEPAPGK" # deleted and re-assigned, continue the loop!


protmtr<-matrix(NA,182,20); for (i in 1:182){ protmtr[i,]<- extractAAC(tg[[i]]);  print(paste(i,"--th round, ends----"))}

# odel 81-st has problem again!
tg[[81]]<-"MSWAKQRVPFLDDDDGEEENDVQDDVDSPVPTRPLVIDEDAEPAAGTSGGLEGGGGDDEDGEDGHALPDLDDDLLLQFEPMLPRVYDLLLPSLDARLNFVNAGQKYAAFLKYVHGDCATCSHGEILREKTQLLTAIVSKLMDINGILEGKDEPAPGK"

tg[[82]]<-"MLRGDSAAKIQERYAELQKRKSHPTSCISTAFTNVATLCRKRYQMMHPELGLAHSCNEAFLPLMAFCGRHRDYNSPEESQRELLFHERLKSALDKLTFRPCSEEQRASYQKLDALTELYRDPQFQQINNFMTDFKKWLDGGFSTAVEGDAKAIRLEPFQKNLLIHVIFFIAVTKIPVLANRVLQYLIHAFQIDFLSQTSIDIFKQKATVFLVPRRHGKTWFIIPIISFLLKHMIGISIGYVAHQKHVSQFVLKEVEFRCRHTFARDYVVENKDNVISIDHRGAKSTALFASCYNTNSIRGQNFHLLLVDEAHFIKKEAFNTILGFLAQNTTKIIFISSTNTTSDSTCFLTRLNNAPFDMLNVVSYVCEEHLHSFTEKGDATACPCYRLHKPTFISLNSQVRKTANMFMPGAFMDEIIGGTNKISQNTVLITDQSREEFDILRYSTLNTNAYDYFGKTLYVYLDPAFTTNRKASGTGVAAVGAYRHQFLIYGLEHFFLRDLSESSEVAIAECAAHMIISVLSLHPYLDELRIAVEGNTNQAAAVRIACLIRQSVQSSTLIRVLFYHTPDQNHIEQPFYLMGRDKALAVEQFISRFNSGYIKASQELVSYTIKLSHDPIEYLLEQIQNLHRVTLAEGTTARYSAKRQNRISDDLIIAVIMATYLCDDIHAIRFRVS" # fixed, deleted the shit space at the end

# haha 153 had problem again

tg[[153]] # outputed "XXXXXDXAPEARQAIRSLTERLYXGGPLTNSKGQNCGYRRCRASGVLTTSCXNTLTCYLKASAACRAAKLQDCT" ------ shit, this fuck shit is really with X amino acids!----just again random AA to it !


tg[[153]] <- "AKDYVDCAPEARQAIRSLTERLYPGGPLTNSKGQNCGYRRCRASGVLTTSCRNTLTCYLKASAACRAAKLQDCT"

# 164 had problem again

tg[[164]]<- "MLRGDSAAKIQERYAELQKRKSHPTSCISTAFTNVATLCRKRYQMMHPELGLAHSCNEAFLPLMAFCGRHRDYNSPEESQRELLFHERLKSALDKLTFRPCSEEQRASYQKLDALTELYRDPQFQQINNFMTDFKKWLDGGFSTAVEGDAKAIRLEPFQKNLLIHVIFFIAVTKIPVLANRVLQYLIHAFQIDFLSQTSIDIFKQKATVFLVPRRHGKTWFIIPIISFLLKHMIGISIGYVAHQKHVSQFVLKEVEFRCRHTFARDYVVENKDNVISIDHRGAKSTALFASCYNTNSIRGQNFHLLLVDEAHFIKKEAFNTILGFLAQNTTKIIFISSTNTTSDSTCFLTRLNNAPFDMLNVVSYVCEEHLHSFTEKGDATACPCYRLHKPTFISLNSQVRKTANMFMPGAFMDEIIGGTNKISQNTVLITDQSREEFDILRYSTLNTNAYDYFGKTLYVYLDPAFTTNRKASGTGVAAVGAYRHQFLIYGLEHFFLRDLSESSEVAIAECAAHMIISVLSLHPYLDELRIAVEGNTNQAAAVRIACLIRQSVQSSTLIRVLFYHTPDQNHIEQPFYLMGRDKALAVEQFISRFNSGYIKASQELVSYTIKLSHDPIEYLLEQIQNLHRVTLAEGTTARYSAKRQNRISDDLIIAVIMATYLCDDIHAIRFRVS" # fixed the shit space 


tg[[166]]<- "MLRGDSAAKIQERYAELQKRKSHPTSCISTAFTNVATLCRKRYQMMHPELGLAHSCNEAFLPLMAFCGRHRDYNSPEESQRELLFHERLKSALDKLTFRPCSEEQRASYQKLDALTELYRDPQFQQINNFMTDFKKWLDGGFSTAVEGDAKAIRLEPFQKNLLIHVIFFIAVTKIPVLANRVLQYLIHAFQIDFLSQTSIDIFKQKATVFLVPRRHGKTWFIIPIISFLLKHMIGISIGYVAHQKHVSQFVLKEVEFRCRHTFARDYVVENKDNVISIDHRGAKSTALFASCYNTNSIRGQNFHLLLVDEAHFIKKEAFNTILGFLAQNTTKIIFISSTNTTSDSTCFLTRLNNAPFDMLNVVSYVCEEHLHSFTEKGDATACPCYRLHKPTFISLNSQVRKTANMFMPGAFMDEIIGGTNKISQNTVLITDQSREEFDILRYSTLNTNAYDYFGKTLYVYLDPAFTTNRKASGTGVAAVGAYRHQFLIYGLEHFFLRDLSESSEVAIAECAAHMIISVLSLHPYLDELRIAVEGNTNQAAAVRIACLIRQSVQSSTLIRVLFYHTPDQNHIEQPFYLMGRDKALAVEQFISRFNSGYIKASQELVSYTIKLSHDPIEYLLEQIQNLHRVTLAEGTTARYSAKRQNRISDDLIIAVIMATYLCDDIHAIRFRVS" # a shit space again




# continue to fuck the loop again? when encountered, manually fuck the spaces?

#again continue
protmtr<-matrix(NA,182,20); for (i in 1:182){ protmtr[i,]<- extractAAC(tg[[i]]);  print(paste(i,"--th round, ends----"))} # cool finally no error !!

# check
head(protmtr) # cool.should be fine.


# I should save the fixed format FASTA file and also the multiple matrice files!!!!-----------both output to local!!!

# (1) export the fixed fasta file
write(tg,data="N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA") # odel error argument !


write(tg,file="N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA") # odel still error shit!

write(tg[],file="N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA") # odel still error shit!

write(tg[[]],file="N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA") # outputed , but shit only saved one line, useless!

writeLines(tg, "N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA") # error !



write.table(tg,"N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA" ) # write.table () worked, but formats are into one-line! format bullshit!

class(tg) # it is a list type!

write.table(tg,
            file = "N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/FixedErrorRstudio182TargetSeqs.FASTA",
            append = F,
            quote = T,
            sep=";",
            eol= "\n"
) # useless shit ! still formats are shits!s



# file Connection another way ---- cool! worked! though with bugs as well

fc<-file("N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/output.fasta"); writeLines(as.character(tg),fc) ; close(fc) # bullshit headers > xxxx seq info were missing, but at least the seqs were output and saparated by lines! -----finally changed name to "FixedErrorRstudio182TargetSeqs.fasta" 



# (2) export the drug 1024 fingerprint matrix combined with label 1/0 column

drugmt<-cbind(lab,drugfp@fpma)
drugmt
drugmt[1,]
head(drugmt[1,])


write(drugmt, "N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/drugFPlabelMatrix.csv") # odel. outputed but very shit strange!!!! many rows and only 4 columns!

# use write.csv directly!
write.csv(drugmt, "N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/drugFPlabelMatrix.csv", sep=",") # messive shit format, but anyway, stored and outputed corrected the data !


# (3) export the 182 protein seqs' amino acids composition matrix 182 rows x 20 columns 
write.csv(protmtr, "N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/182proteinProtRAACMatrix.csv") # same above, messive thoughs

# should try to re-load it into the R studio, in case that like the shit fasta file, fixed in the Rstudio, save, but failed to reload it normally as fasta variable in Rstudio

reld<-read.csv("F:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/182proteinProtRAACMatrix.csv")
# check if reloaded normally
class(reld) # data.frame type, not matrix, but does not matter !
dim(reld) # 182 x 21! why not 20? with label?
tail(reld) # checked, the 1st column was the ID/ order/ index !



#  combination into full  matrix for machine learning !
cbind(drugmt,protmtr) # cool! no error!
dim(cbind(drugmt,protmtr)) # cool! 182 x 1045

fulmat<-cbind(drugmt,protmtr)
dim(fulmat) # cool!

# also export it to local file 
write.csv(fulmat, "N:/1 prjsRes/开题强制hpvDrugRepo呵呵6jul2020/fullMatrix182PairsDrugTarget4ML_AACdescFeat.csv")


# random forest !  go!
install.packages("RandomForest") # bullshit uppercase R was wrong!
install.packages("randomForest") # cool

library(randomForest) # cool called
rfNews() # updates for new version

randomForest(x = fulmat[1:146,2:1045], 
             y = fulmat[1:146,1],  
             xtest=fulmat[147:182,2:1045], 
             ytest=fulmat[147:182,1], 
             ntree=500 ) # odel outputed, but with warning : The response has five or fewer unique values.  Are you sure you want to do regression?????? regressin? ------fuck shit not showing performance? 



randomForest(x = fulmat[1:146,2:1045], 
             y = as.factor(fulmat[1:146,1]),  
             xtest=fulmat[147:182,2:1045], 
             ytest=as.factor(fulmat[147:182,1]), 
             ntree=500 ) # added as.factor, cool outputed differently! as below:

"""
Call:
 randomForest(x = fulmat[1:146, 2:1045], y = as.factor(fulmat[1:146,      1]), xtest = fulmat[147:182, 2:1045], ytest = as.factor(fulmat[147:182,      1]), ntree = 500) 
Type of random forest: classification
Number of trees: 500
No. of variables tried at each split: 32

OOB estimate of  error rate: 13.7%
Confusion matrix:
0  1 class.error
0 67  6  0.08219178
1 14 59  0.19178082
Test set error rate: 41.67%
Confusion matrix:
  0  1   class.error
0 4 14  0.77777778
1 1 17  0.05555556  =====odel totally 36 validation set? why? not 34?? each 17? ---------shit, i am wrong , each (pos and nega) validation set was 18, not 17. totally was 36 not 34 ! I was wrong!




above meaning that 00 is the true negative -- 4 
11 is the true positive ! 17

usually column names are horizontally written ! so the horitonzal header 0 1 are indicating that predicted-NO and predicted-YES, respectively. 

so the vertical 
0
1 are actually Y/N, 1/0 , T/F . 

 so  , 14 here is the value of false positive, as they were predicted to be true but actually they are not!

so here,  true postive rate= true positive instances / all positive cases = 17/18 = 0.94444

true negative rate= number of true negative / all negatives = 4 /18 = 0.2222

false negative rate = 1/18 = 0.05555

false positive rate= 11/18 = 0.61111

"""



# trial without testing set ? will return predictive performance rate?
randomForest(x = fulmat[,2:1045], y = fulmat[,1], ntree=500 ) # odel shit, went back to regression again!

randomForest(x = fulmat[,2:1045], y = as.factor(fulmat[,1]), ntree=500 ) #  heiyo !!!! odel !!! outputed below shits!! ----------- so  , randomForest itself is self-attached for k-fold cross validation??? or??? do I need to use tag_label_shuffle to generate multiple models like previous BRCA1-ML work???

"""
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

calculation:
true positive rate: 75/total positive(182/2=91) = 75/91 = 0.82

true negative rate: 81/91 = 0.88

false positive rate = 10/91=  0.1098901 ~ 11%

false negative rate = 16/91= 0.1758242 ~ 17.58%
"""

# plot above two set of confusion matrices!?
install.packages("ROCR") 
library(ROCR) # cool , installed, loaded!


performance(prediction(c(0,0,1,1),c(1,0,1,1)), "auc") # heiyo. no error ! can directly plot?

plot(performance(prediction(c(0,0,1,1),c(1,0,1,1)), "auc")) # odel cannot plot directly

ac<-performance(prediction(c(0,0,1,1), c(1,0,1,1)), "auc"); plot(ac) # shit still error!



# create two var storing labels.

labl<-list(); pred<-list()
# firstly, ( i was stupid, ) the actual label set/validation labels' number were known , i should firstly assign all the known 0 and 1 values to the known label var !

# secondly,  true positive 17, total positive 18. So , to add 17 to "1" to the "pred" var and then 18 "1" to labl var

pred<-rep.int(1,17); print(pred); length(pred) # cool!
labl<-rep.int(1,18); 

# real labels are known in advance, so can set it first, the latter 18 elements of the real labels are 0, as i pre-setted a balanced validation set, so
labl<-rep.int(0,18); tail(labl); length(labl) # cool

# the 18th element of pred must be "0"! so
pred[18]<-0; print(pred[18]) # cool! it was the only element of false negative

length(labl) # 18 only . hehe add to full 36.
labl[19:36]<-0; print(labl); length(labl);  # cool 36 elements.
labl[17:19]

# and then fuck the 14 false positive results.
length(pred) # 18
pred[19:32]<-1; pred[33:36]<-0; print(pred)


# notice that in processing of 17 true positive there was already 0(predicted) and 1(actual, i.e.,false negative) at 18th element --- yep, the only false negative one

library(ROCR) # switched into another PC machine's Rstudio. hehe shit, so has to recall library again

pd<- prediction(pred, labl);   pdd<-performance(pd, "auc")  ; pdd@y.values # odel 0.58333!! only sucks!

ptd<-performance(pd,"tpr","fpr"); plot(ptd)# plotted, but result sucks!

plot(ptd, abline=c(0,1),add=T) # shit, no diagonal was added!
plot(ptd, a=0.5 ,add=T) # shit, nothing was added!
plot(ptd, x=0.5 ,add=T) # shit, nothing was added!
plot(0.5 ,add=T) # shit, error
plot(x=0.5 ,add=T) # shit, error
plot(x=0.3, y=0.4 ,add=T) # shit, just a point  , the curve disappeared!

# plot 2nd model's 
# prepare the labelset first: 182 in total, half positive half negative
l2<-numeric()
l2[1:91]<-1; l2[92:182]<-0; length(l2); l2[91:93] # cool! has to create the var<-numeric();

# and then the predictive results'  :  
# first:  75 true positive -- cannot write other because the order of actual label set was determined, and predictive set has to match the pre-setted actual label set --- predictive's 1s should match actual alabel set's 1s
p2<-numeric(); p2[1:75]<-1; 
# second: 16 false negative
p2[76:91] <- 0; length(p2) # cool 91. hehe, 
p2[75:91] # cool, correctly output
# third : 10 false positive and  fourth : 81 true negative:
p2[92:101]<-1; p2[102:182]<-0; print(p2[91:182])

pn<-prediction(p2,l2)

print( performance(pn, "auc")@y.values) # odel not bad ! 0.8571489!!!

pk<-performance(pn, "tpr","fpr"); plot(pk)   

"""  
2nd random forest model, without self labelled testing/validation set's resultant confusion matrix
   0  1
0 81 10 
1 16 75


1st random forest model, self labelled testing/validation set's resultant confusion matrix

Confusion matrix:
  0  1                   predicted false    predicted true
0 4 14    actually false        4                 14
1 1 17    actually true         1                17 

"""



### odel DO I still need to export the random forest model(s)? into local file ????????? shit, it was in the HP omen 17, not in other PC!! need to retrain and load many data if want to get the model here !
###  or I should save R.data and make it portable ?? store in USB disk or mobile harddisk?
### 

# continue try other further protein feature descriptors for randomForest loading training ??? 


