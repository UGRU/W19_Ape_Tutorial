####TUTORIAL: INTRODUCTION TO THE ANALYSIS OF DNA SEQUENCE DATA FROM BOLD USING R----

#last updated March 8, 2019 by Sally Adamowicz (sadamowi@uoguelph.ca)

#Acknowledgements: Some lines below were inspired by tips or work by PhD candidates Jacqueline May and Matthew Orton. As well, thank you to the UGRU group for many helpful tutorials, resulting in my near-complete conversion to the tidyverse.

####PART 1 - OVERVIEW: WHAT IS COVERED IN THIS TUTORIAL?----

#Introduction: The Barcode of Life Data Systems (BOLD) database (www.boldsystems.org) houses DNA sequence data from standardized gene regions. These data are commonly used for the identification of unknown specimens as well as for biodiversity discovery and exploration. One advantage of BOLD over NCBI is the structure of the metadata. Many records on BOLD have GPS coordinates in a standardized format, which facilitates geographic analysis. When combined with molecular taxonomic units to overcome the taxonomic impediment, the geo-referenced specimen and sequence data provide a powerful resource for new research on the phylogeography and macroecology of understudied organisms.

#This tutorial integrates several core biological analysis objectives together with several useful functions and programming tools.

#Biological learning outcomes: By the end of this tutorial, we should be able to: download data from BOLD directly into R using the API tool, apply filters to the data, explore nucleotide composition, perform DNA sequence alignment, cluster our data into Molecular Operational Taxonomic Units (MOTU), explore how our choice of clustering threshold impacts the number of MOTU generated from a set of sequences, and build a neighbour-joining tree using DNA barcode sequence data.

#Programming learning outcomes: By the end of this tutorial, we should be able to: filter our data in a logical way using piping and helpful functions from the tidyverse; access data from dataframes, lists, and other list-type data objects in R using indexing; perform simple examples of string matching using regular expressions; use functions from a few packages from the Bioconductor repository; and explore the use of the "for loop" vs. "while loop" in the context of iterating over clustering thresholds.

#Note about iteration: Generally, as we have seen in prior UGRU tutorials, the "apply" family of functions is preferred if we are performing exactly the same operation across the elements of a data object. For example, we would use the function lapply() for to iterate over the elements of a list. However, programming loops ourselves is also useful when we want to change something with each iteration. So, below we explore the for loop and while loop for iterating across various clustering thresholds. I would also be curious if you have a better solution to the below problems (e.g. perhaps using the map functions we learned about last week?), and I'd be glad to see your solutions. Please do share.

#Feedback: I welcome your feedback on this tutorial. It is my goal for the script to be accessible for all user levels, and I would appreciate any comments that you have. 


####PART 2 - LOAD PACKAGES----

#Version: This script was prepared and run using R Version 3.5.2 (Eggshell Igloo) and RStudio Version 1.1.423. Most functions below should work fine on other recent versions of R (but you may wish to update if you want to explore the dendextend package).

#TIPS in case you have trouble with installing packages: If you haven't already installed these, you would need to install the packages before you can load them. You would remove "#" and would need to run the installation lines. Installation only needs to be performed once, unless you wish to reinstall R or update a package. By contrast, you would need to load the libraries each time you start a new session. As well, I suggest to read error messages carefully. For example, there could be dependencies you need to install. It is also important to be sure you are seeking each package from the correct repository, e.g. CRAN or Bioconductor or another source. Finally, you can try changing the CRAN mirror site being used if you are having difficulties with a download and installation for a package that ought be on CRAN. You would do that using the following function. You would uncomment this and then, in the console, select the CRAN mirror you want to try. A final potential solution if you are having trouble with package installations is to reinstall R and RStudio. Generally, it can be helpful to use updated versions in order to be able to use cutting-edge packages and the newest versions of commonly-used packages.

#Uncomment this line if you want to choose a different CRAN mirror site.
#chooseCRANmirror()

#Packages from the CRAN repository that we will use for this tutorial

#We will use functions from tidyverse packages for filtering our data and searching for strings. We will also make use of the pipe from the maggitr package of the tidyverse suite.
#install.packages("tidyverse")
library(tidyverse)

#ape is a core package for working with phylogenies in R
#install.packages("ape")
library(ape)

#This package provides additional functionality for visualizations of dendrograms. 
#install.packages(dendextend)
#library(dendextend)

#Packages from Bioconductor, which is a repository that focuses on curated packages for bioinformatics. You would need to uncomment and run all of these lines if you don't have these packages yet. Bioconductor is a very helpful resource, housing >1500 packages.

#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("muscle")
#biocLite("msa")
#biocLite("DECIPHER")

library(Biostrings)
library(muscle)
library(msa)
library(DECIPHER)

#Note about the Bioconductor packages: Biostrings is a very useful package for manipulating biological sequence data and for converting data objects between classes. It also has functionality for exploring sequence features (e.g. nucleotide composition, k-mer frequencies) and for analysis (e.g. mapping sequence reads to a reference genome). We will use muscle here for DNA sequence alignment. The package msa has additional functionality for building and viewing alignments. DECIPHER is one of my favourite packages; beyond this tutorial, I would recommend to check out its diverse functionality, which spans clustering, alignment, primer design, among other topics.

#Vignettes. Bioconductor packages typically have a vignette, which provides a brief overview of some of the most important and typical features of the package.

#For example, we can see what vignettes are available for the Biostrings package as follows:
vignette(package = "Biostrings")

#We can then bring up the selected vignettes we are interested in, for example:
vignette("BiostringsQuickOverview", package="Biostrings")
vignette("MultipleAlignments", package="Biostrings")

####PART 3 - IMPORT AND CHECK CONTENTS OF DATASET----

#I chose the genus Cotesia for this tutorial for a couple of reasons: it is a genus of parasitoid wasp with interesting metadata you may wish to explore (e.g. the host organism is available for some of the specimens), and there is more than one DNA marker type available. We can explore differences in a mitochondrial vs. a nuclear marker. Also, much of the data on BOLD for this genus is generated by Alex Smith and colleagues, and so I thought this taxon would be of interest to department members. I recommend this interesting paper that includes this taxon: Smith et al. 2008. Extreme diversity of tropical parasitoid wasps exposed by iterative integration of natural history, DNA barcoding, morphology, and collections. PNAS. 105(34): 12359-12364. doi: www.pnas.orgcgidoi10.1073pnas.0805319105

#You can choose to set your working directory for the session in RStudio. e.g. Session -> Set Working Directory -> To Source File Location (which would set your working directory to the same folder from where you opened this script file.)

#Dataset. This tutorial uses a dataset which was downloaded from BOlD. I suggest to use the file provided, so that we are using exactly the same version and get the same results. However, you may also choose to download this file from BOLD directly for yourself to explore an updated dataset and to see how the download from the API works. The data file was obtained September 19, 2018 using the following line of code:

#Cotesia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cotesia&format=tsv")

#Note that we can jointly call data based upon both taxonomic and geographic criteria if we want to. The above line asks for data only based upon taxonomy. You can read more about BOLD's API tools here: http://www.boldsystems.org/index.php/resources/api?type=webservices. I suggest to try downloading and playing with datasets that interest you.

#We can next write the file to hard disk, as it is helpful to save a copy of the original data prior to manipulation and analysis so we can always go back to the original.
#write_tsv(Cotesia, "Cotesia_BOLD_data.tsv")

#Import Cotesia_BOLD_data.tsv file. Note that we are using the function read_tsv (from the tidyverse package readr). The readr functions tend to have better data parsing functionality than the base R functions. You can read more in "R for Data Science" by Wickham and Grolemund, available at https://r4ds.had.co.nz/
Cotesia <- read_tsv("Cotesia_BOLD_data.tsv")

#Checking class and viewing data summary on the screen. 
class(Cotesia)
summary(Cotesia)

#We can see we have a tibble (a tidyverse style of dataframe, see: https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html). So, we have a rectangular data object, which has various types of data in the columns. A nice thing about BOLD is that the data are available for download in quite tidy format. The individual specimen records are in rows. Variables are in columns. Data values are in the cells. (We can read about tidy data here: https://r4ds.had.co.nz/tidy-data.html). Specimen records from BOLD may or may not have an associated DNA sequence.

#Checking variable names and we can view the data in RStudio by clicking on the name Cotesia in the environment window.
names(Cotesia)

####PART 4 - FILTERING AND FORMATTING DATA----

#Seeing what markers are available in this dataset
unique(Cotesia$markercode)

#Next, let's obtain count of the sequences by markercode. Our end goal is: a simple dataframe with a list of the markercodes, the counts of the numbers of records having each markercode, arranged in descending order by count. Basically, we want to know: what genes have a decent sample size for an analysis for this taxonomic group?

#We can do this easily using tidyverse functions. We first specify that we will assign the results to count.by.marker. After the assignment operator, we specify we are working with the Cotesia tibble. When using piping, we can thereafter use only the column names and not have to refer to the dataframe name over and over (so we don't have to use the dataframe name and $ over and over again. With piping, the output from one step serves as the input for the next step. After specifying the dataframe we are working with, we group the rows by markercode. i.e. We are grouping the records together by gene region. We then can use the length() function to count the number of records having each markercode. Processid is a BOLD identification code for each sample. We declare n as those lengths (i.e. counts by markercode). Note that we are using "=" within a function call, and so n will not exist in our workspace after we run this. Passing the counts as an argument to summarize() will create a new dataframe for us. If we were to use mutate() a new column would be added to our existing dataframe. Here, we want a simple dataframe of just the the markers and record counts. Then, we arrange the new dataframe in descending order by the counts and add a print statement so that we can see the results immediately on the screen. Piping allows us to have streamlined code and a very logical flow. We also avoid creating unneeded interim objects. However, it is important to check that the code is doing what we want.

count.by.marker <- Cotesia %>%
  group_by(markercode) %>%
  summarize(n = length(processid)) %>%
  arrange(desc(n)) %>%
  print()

#We can see that most of the data are COI-5P (the 5' end of the cytochrome c oxidase subunit I mitochondrial gene). There are 96 records missing the specification of the markercode. There are 82 records having markercode 28s (the nuclear large ribosomal RNA subunit). The other genes have a small sample size. So, we will consider COI-5P and 28S to be candidate markers for our analysis.

#Next, we will filter the dataset to retain those records having a COI-5P markercode. Note we need to use "==" for exactly equal. We also want to filter out records lacking a nucleotide sequence. Why might a record lack a sequence? Well, the sequencing may not have been performed yet. Or, a particular gene region was not attempted. Or, the sequencing failed. Or, a sequence was obtained but then subsequently deleted due to the discovery of a contamination event. So, there are many specimen records on BOLD that lack sequence data. Here, we are doing this using a "regular expression". A regular expression is one or a series of characters that will allow us to define and search for a precise pattern in character data (i.e. string data). We are looking to see if there is an A, C, G, or T in each cell in the the nucleotides column. The square brackets specify we are matching A or C or G or T. If we instead wanted to searh for the pattern "ACGT" as a sequence, we would omit the square brackets. We are assigning our COI dataset to a new tibble called Cotesia.COI.

Cotesia.COI <- Cotesia %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[ACGT]"))

#Checking that we have only COI sequences remaining in Cotesia.COI
unique(Cotesia.COI$markercode)

#Similarly, we are next creating a tibble containing records having 28S sequences.

Cotesia.28S <- Cotesia %>%
  filter(markercode == "28S") %>%
  filter(str_detect(nucleotides, "[ACGT]"))

#Recommended reading: I recommend to read the chapter in "R for Data Science" by Wickham and Grolemund on regular expressions (Chapter 11). That is a very helpful, accessible chapter that provides a great introduction to regular expressions. Regular expressions are very helpful for working with string data, such as DNA, RNA, or protein sequences.

#Also, I enjoyed reading this humorously-titled yet also helpful post on regular expressions: http://web.archive.org/web/20090226052234/http://immike.net/blog/2007/04/06/the-absolute-bare-minimum-every-programmer-should-know-about-regular-expressions/
#The basics for regular expressions are similar across programming languages, but we can look up patterns in our language of choice to be sure. Also, as always, we should check that code is doing what we want by checking the performance on test cases.

#Now we can remove our original dataframe Cotesia to clean up our workspace, as we aren't planning to work with the other markers. Note we can easily recreate Cotesia by running the above read_tsv command. We did NOT delete the original tsv from hard disk by doing this. We are just removing Cotesia from our R environment.

rm(Cotesia)
rm(count.by.marker)

#Depending upon the nature of our study, we would of course need to consider applying different filters to the raw data download.

#Next, we will change the format of the nucleotide data to be suitable for downstream analysis. Here, we will change the sequence data to a DNAStringSet object class, using a function from the Biostrings package.
Cotesia.28S$nucleotides <- DNAStringSet(Cotesia.28S$nucleotides)
Cotesia.COI$nucleotides <- DNAStringSet(Cotesia.COI$nucleotides)

#So, how did we know we needed to do that data conversion step? For steps we want to do downstream, we can check what type of data we need for functions we want to use. For example, we can bring up the documentation for the muscle() function from the muscle package as follows:
?muscle::muscle
#We can then see that this function takes a stringset object. It is always essential to know what kind of data class we have and what kind of data class(es) a particular function will accept as arguments. In our case, we have DNA sequence data, and so we converted to a DNAStringSet object using the DNAStringSet() function.

#Further readings about data structures. Packages available through the Bioconductor system use the S4 system of object-oriented programming.
#http://adv-r.had.co.nz/OO-essentials.html
#http://adv-r.had.co.nz/S4.html
#https://master.bioconductor.org/help/course-materials/2017/Zurich/S4-classes-and-methods.html
#https://bioconductor.org/help/course-materials/2013/CSAMA2013/friday/afternoon/S4-tutorial.pdf

#What does this mean from a practical point of view? The StringSet objects are similar to lists. Therefore, we would use double square brackets [[]] to retrieve an element for analysis; i.e. if  we want to work with the data from an element in a list object. By contrast, for subsetting into another list object, we would use a single set of square brackets []. We will use these indexing principles for working with our sequence data in downstream steps.

#So, we now have two tibbles filtered from our original data download, one for COI Sequences and one for 28S sequences, and are reading for some further data analysis.


####PART 5 - CALCULATING NUCLEOTIDE COMPOSITION AND USING GGPLOT----

#Checking the data type of the nucleotides
class(Cotesia.COI$nucleotides)
class(Cotesia.28S$nucleotides)

#Next, we will use some functions from the Biostrings package to explore the sequence data.

#For accessing the help for Biostrings function, we leave off the ()
?letterFrequency

#This will calculate the nucleotide frequency of the specified letter, in this case A. Here, we want to look at the results and so we are coercing the nucleotide frequencies to a data frame object and assigning the results to a named dataframe.
Cotesia.COI.AFreq <- as.data.frame(letterFrequency(Cotesia.COI$nucleotides, letters = "A"))

#We can now easily have a look in the viewer by clicking on the new dataframe name in the Global Environment viewer (upper right corner in RStudio).

#This next line will calculate the nucleotide frequency of the specified letters, A or T. So, this calculates the frequency at which an A or T is observed across positions. Run this and have a look.
Cotesia.COI.ATFreq <- as.data.frame(letterFrequency(Cotesia.COI$nucleotides, letters = "AT"))

#By contrast, the next line calculates the frequency of each nucleotide separately.
Cotesia.COI.NucFreq <- as.data.frame(letterFrequency(Cotesia.COI$nucleotides, letters = c("A", "C", "G", "T")))

#Note the default for the function letterFrequency is to return counts. The counts will depend upon both nucleotide composition and sequence length.

#We could calculate the AT frequency in a different format that is easier to interpret, as a proportion of the total sequence length. Here, we are using mutate to add a column onto the end of Cotesia.COI.NucFreq. We don't want to repeat the dataframe name numerous times. Here, we are using piping to specify the source dataframe once at the beginning. Alternatively, the dataframe name can be specified as the first argument of mutate.
Cotesia.COI.NucFreq <- Cotesia.COI.NucFreq %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))

#Compare with the following:
Cotesia.COI.NucFreq1 <-  mutate(Cotesia.COI.NucFreq, ATproportion = ((A + T) / (A + T + G + C)))

#Are these the same?
all.equal(Cotesia.COI.NucFreq, Cotesia.COI.NucFreq1)

#We can append the processid's so that we can link back the nucleotide frequencies to specific variables if we want to. We want to do this right away before working further with the data.
Cotesia.COI.NucFreq$processid <- Cotesia.COI$processid

#Coding challenge! Make a dataframe like Cotesia.COI.NucFreq and append the processid's in one step.

#This one was a little tricky, and so don't worry if you didn't get this one. Here is one potential solution. Do check it. With piping, the output from one step is the input to the next step. Typically, with piping, we would be specifying a single dataframe after the assignment operator and working with that. Whichever method you choose, it is important to check that it worked as expected. As well, it would be important to append the unique identifiers right away and not risk the data getting dissociated.

Cotesia.COI.NucFreq2 <- as.data.frame(letterFrequency(Cotesia.COI$nucleotides, letters = c("A", "C", "G", "T"))) %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C))) %>%
  bind_cols(processid = Cotesia.COI$processid)

#NOTE: bind_cols() is the tidyverse function simlar to cbind()

#Now, let's look at the histogram of AT proportions.
hist(Cotesia.COI.NucFreq$ATproportion)

#Mean AT proportion
mean(Cotesia.COI.NucFreq$ATproportion)

#Let's do the same for 28S. Doing counts of letters.
Cotesia.28S.NucFreq <- as.data.frame(letterFrequency(Cotesia.28S$nucleotides, letters = c("A", "C", "G", "T")))

#Calculating AT proportions.
Cotesia.28S.NucFreq <- Cotesia.28S.NucFreq %>%
  mutate(ATproportion = ((A + T) / (A + T + G + C)))

#Note we can append the processid's so that we can link back the nucleotide frequencies to specific samples if we want to.
Cotesia.28S.NucFreq$processid <- Cotesia.28S$processid

#viewing histogram. (We will come back to the unusual low value in a minute.)
hist(Cotesia.28S.NucFreq$ATproportion)

#calculating mean
mean(Cotesia.28S.NucFreq$ATproportion)

#Are the mean AT proportions different between COI and 28S in the genus Cotesia? We can perform a t-test.
t.test(Cotesia.28S.NucFreq$ATproportion, Cotesia.COI.NucFreq$ATproportion)

#Is the result still the same when we omit an apparent outlier? (We will address that specific value further below.)
t.test(Cotesia.28S.NucFreq$ATproportion[-c(65)], Cotesia.COI.NucFreq$ATproportion)

#Code to prepare histogram figures using ggplot. Note that above we used hist(), which works fine for simple plots and can be somewhat customized (e.g. changing axis labels, colours, etc.). ggplot is much more flexible for creating graphics. I suggest to read Chapters 1 and 5 of "R for Data Science" by Wickham and Grolemund.

#Just COI. Here, I used the defaults for bin size, etc. Note that we don't have to use the names of all the arguments if we use your arguments in order. e.g. You don't have to type "data = " and "mapping = ".
ggplot(data = Cotesia.COI.NucFreq, mapping = aes(Cotesia.COI.NucFreq$ATproportion)) + 
  geom_histogram() +
  labs(title="Histogram of Cotesia COI Nucleotide Composition", x="AT Composition", y="Count")

#Plotting COI again, but see how we are shifting bin size and the x-axis limits. In increments, we are setting up the plot also to receive the 28S data, so that we can see what the different additions to the code do.
ggplot(data = Cotesia.COI.NucFreq, aes(Cotesia.COI.NucFreq$ATproportion)) + 
  geom_histogram(breaks=seq(0.3, 0.8, by = 0.01), 
                 col="black", 
                 fill="red", 
                 alpha = .2) +
  labs(title="Histogram of Cotesia COI Nucleotide Composition", x="AT Composition", y="Count") +
  xlim(c(0.3,0.8)) +
  ylim(c(0, 2000))

#Joining together the COI and 28S NucFreq dataframes. First, we want to add the markercode as a new column to each NucFreq dataframe so that there will be a clear indicator of the markercode when we do the data join. Then, we can use bind_rows() to add the rows from the 28S nucleotide frequency file onto the row of the COI nucleotide frequency file. Remember that we can use the "Data Wrangling with dplyr and tidyr" Cheat Sheet to look up which function we want for a particular type of data manipulation. Here, bins_rows() does the trick.
Cotesia.28S.NucFreq$markercode <- "28S"
Cotesia.COI.NucFreq$markercode <- "COI"
Cotesia.COI.28S.NucFreq <- bind_rows(Cotesia.COI.NucFreq, Cotesia.28S.NucFreq)

#Now, ready to plot histograms on the same plot, coloured by markercode, allowing us to compare visually the distribution of AT content between these marker codes.
ggplot(data = Cotesia.COI.28S.NucFreq, aes(Cotesia.COI.28S.NucFreq$ATproportion)) + 
  geom_histogram(mapping = aes(colour = markercode, fill = markercode),
                 breaks=seq(0.3, 0.8, by = 0.01), 
                 alpha = .3) +
  labs(title="Histogram of Cotesia COI Nucleotide Composition", x="AT Composition", y="Count") +
  xlim(c(0.3,0.8)) +
  ylim(c(0, 2000))

#You can try playing around with the different parameters to see what happens. For example, we can vary the setting for "alpha" to create different levels of transparency. Some transparency is helpful for overlaying histograms for different data on top of one another. Here, we will try plotting the histograms for both of our markers we are analyzing on the same plot, to enable comparison.

#Note that there are very different sample sizes between markers. How should we deal with this to make the comparison between COI and 28S easier to see?

#Coding challenge!! Change the type of plot from a histogram to a density. Hint, we can use geom_freqpoly() rather than geo_histogram(). Also, you would need to change the y-axis label and range. If you make a density plot, be sure to record the sample size information somewhere (such as in the figure title or figure legend) so that the sample size information isn't lost. Sample size information is important for data interpretation.

#Potential answer:
ggplot(data = Cotesia.COI.28S.NucFreq, mapping = aes(x = Cotesia.COI.28S.NucFreq$ATproportion, y = ..density..)) + 
  geom_freqpoly(mapping = aes(colour = markercode, fill = markercode),
                breaks=seq(0.3, 0.8, by = 0.01),
                alpha = 1,
                size = 1.5) +
  labs(title="Density Plot of Cotesia Nucleotide Composition", x="AT Composition", y="Density") +
  xlim(c(0.3,0.8)) +
  ylim(c(0, 70))

#Note you can play with the colours, line width, etc. to optimize the plot way you want. This plot more clearly shows the different distributions in AT composition, compared to the histogram. The density under each curve is the same.

#Moving on to considering a more complex set of sequence features: k-mers. Here, we are calculating the counts of short oligonucleotide sequences of length 4. Counts of all possible 4-nucleotide strings are tabulated. k-mers can be a very useful sequence feature. For example, k-mers are being used for sequence classification in metagenome studies. k-mers are also used during sequence alignment by some algorithms (see paper by Edgar 2004 on MUSCLE algorithm). We could choose a different length for our k-mers. Here, we are choosing 4. You could choose a different number and see what happens.
oligo <- oligonucleotideFrequency(Cotesia.COI$nucleotides, 4)

#But... let's return to that unusual 28S value...

#What is the minimum value of ATproportion in 28S?
min(Cotesia.28S.NucFreq$ATproportion)

#Which sample exhibits the minimum value? this will give row number.
which.min(Cotesia.28S.NucFreq$ATproportion)

#Or, we could have used indexing to obtain the processid (unique identifier for each sample, used by BOLD) of the sample with the minimum AT proportion.
Cotesia.28S$processid[which.min(Cotesia.28S.NucFreq$ATproportion)]

#Checking that there was just one very low value. What is the min value after leaving out row 65?
min(Cotesia.28S.NucFreq$ATproportion[-65])

#Converting that specific sequence to a character vector
to_blast <- as.character(Cotesia.28S$nucleotides[65])

#checking class
class(to_blast)

#printing to screen
to_blast

#We can use this sequence to explore the BLAST tool offered through NCBI. This sequence is possibly a misidentification. This sequence BLASTEd mainly to other members of order Hymenoptera. This sequence could be dropped or explored further prior to further analysis, such as phylogenetic analysis, depending upon the goals of a particular study. In such a case, I would definitely recommend to code in data exclusion and document your reasons for this in the commenting of our script (and into the manuscript, if suitable). As well, for phylogenetic analysis for the genus Cotesia using 28S, I suggest to perform the multiple sequence alignment (MSA) again after omitting the outlying sequence. A sequence that is very different from the target dataset can influence the MSA.

#We can objects we don't need further from our environment using rm(). For example:

rm(Cotesia.COI.AFreq, Cotesia.COI.ATFreq, Cotesia.28S.NucFreq, Cotesia.COI.28S.NucFreq, Cotesia.COI.NucFreq, Cotesia.COI.NucFreq1, Cotesia.COI.NucFreq2, oligo)

#Conclusion from this part: Using calculated values, a statistical test, and graphical exploration, we have explored nucleotide composition between two markers in the parasitoid wasp genus Cotesia. We have seen that composition varies between 28S and COI. The mitochondrial marker is actually highly AT biased, which is a common pattern. We would wish to consider this in our choice of a genetic distance metric for phylogenetic analysis.

####PART 6 - ALIGNMENT----

#Next, we will use the muscle algorithm from the muscle package to perform a multiple sequence alignment. While this is one of the most commonly-used alignment algorithms, there are many alignment algorithms in existence to consider for your specific study. The original algorithm is described in: Edgar 2004.  MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research, 2004, Vol. 32, No. 5, 1792-1797. DOI: 10.1093/nar/gkh340

#For now, we will work with the 28S dataset. You may wish to explore the COI dataset on your own. Depending upon your study goals, you may not want to wait for an alignment that takes a long time to run. One might, for example, randomly select one sequence per BIN (Barcode Index Number) or species or MOTU (Molecular Operational Taxonomic Unit) for alignment and phylogenetic analysis.

#If we want a fast preliminary alignment, we can set the argument maxiters (maximum iterations) to a small value, such as 2. See Edgar paper and also the online muscle documentation: http://www.drive5.com/muscle/muscle_userguide3.8.html. The output of the function will be an object of type MultipleAlignment.
Cotesia.28S.alignment <- muscle::muscle(Cotesia.28S$nucleotides, maxiters = 2)

#The output will be of class DNAMultipleAlignment (which we can later convert to other formats, such as DNAStringSet, as needed for downstream analysis).
class(Cotesia.28S.alignment)

#Assigning sequence names in the alignment to be the processid. The processid is a unique identifier assigned by BOLD.
rownames(Cotesia.28S.alignment) <- Cotesia.28S$processid

#Outputting the alignment to FASTA format, so that we can readily view the alignment (e.g. in a software such as MEGA). There are various ways to do this. To use the function writeXStringSet, note that we had to convert our alignment to a DNAStringSet object, a conversion which is nested in the below line rather than creating a new object in R at this time.
writeXStringSet(DNAStringSet(Cotesia.28S.alignment), file = "msaCotesia28S.fas", format = "fasta")

#The package msa offers interesting functionality for viewing alignments, but requires installation of LaTeX. I suggest to check out: https://rdrr.io/bioc/msa/f/inst/doc/msa.pdf. For now, I suggest to use the user-friendly software MEGA: https://www.megasoftware.net/.

#When we check out the alignment, we see that the sequence in position 65 (GBMIN142670-18) is quite different from others. When BLASTing this, we noticed that most hits are to a different genus. Could this be a case of misidentification? We will check out the genetic distances further below using a distance matrix and dendrogram.

#Note that there is only partial help directly available through R regarding the implementation of the muscle aligment algorithm. The following doesn't bring up all the options.
?muscle::muscle

#However, the help does link to the MUSCLE website, which explains the options in more detail. For example, for the maxiters option: "You can control the number of iterations that MUSCLE does by specifying the -maxiters option. If you specify 1, 2 or 3, then this is exactly the number of iterations that will be performed. If the value is greater than 3, then muscle will continue up to the maximum you specify or until convergence is reached, which ever happens sooner. The default is 16. If you have a large number of sequences, refinement may be rather slow." From: http://www.drive5.com/muscle/muscle_userguide3.8.html Note that the default mentioned in this quote is for the online muscle tools. Below, we will also want to have a look at the default options used specifically as implemented in the R package muscle.

#removing preliminary alignment
rm(Cotesia.28S.alignment)

#Let's try out different analysis settings. For example, what happens when we dramatically change the gap penalties?

#First, the below command does not alter the default for maxiters to allow the alignment to improve over multiple iterations if needed. Also, we are leaving the gap opening penalty and other parameters at the default values as well. We are requesting a log file so that we can see the defaults used in this R implementation of the muscle algorithm. We can open the log file, which was saved in the working directory. By looking at that file, we can see that the default gap opening penalty is -400 and the default for maxiters is 8. This step is helpful because we want to understand the specifics of how the algorithm is implemented in this package.
Cotesia.28S.alignment1 <- DNAStringSet(muscle::muscle(Cotesia.28S$nucleotides, log = "log.tx", verbose = T))

#Here, we are assessing the length of the alignment, using the first element. We will use indexing with double square brackets here, similar to the case of lists, to dive in to retrieve the data from the first element to be able to analyze it. If we used a single set of square brackets, that would be for subsetting and would return another list object, rather than pulling out the data for us to work with. So, here we are seeing how many alignment positions are in the first sequence in the alignment.
length(Cotesia.28S.alignment1[[1]])

#The alignment is 591 positions long.

#Now, let's try a different alignment, with a very different gap opening penalty. This is penalizing the opening of gaps more heavily for the calculation of the alignment score.
Cotesia.28S.alignment2 <- DNAStringSet(muscle::muscle(Cotesia.28S$nucleotides, gapopen = -10000))

#Checking length of the alignment
length(Cotesia.28S.alignment2[[1]])

#Here, the length is much shorter, 533 positions, as gaps are penalized much more heavily. Basically, our choice of a heavy gap opening penalty is forcing sequence sections together that are more disparate, in comparison with alignment1.

#So, we can see that the total alignment lengths from the same starting sequences can be different depending upon the gap penalty. We could also count the number of gaps.

#Counting the number of gaps in sequence 1 of alignment 1. We are matching to a specific character, "-".
str_count(Cotesia.28S.alignment1[[1]], "-")

#Now let's look at the sequence in row 65, which was an unusual sequence according to the above analysis looking at nucleotide composition. Also, we had blasted the sequence and concluded that it could have been a misidentified specimen actually belonging to a different genus.
str_count(Cotesia.28S.alignment1[[65]], "-")

#Now, we will look at these metrics for these specific sequences in alignment 2. The first sequence has just one gap in alignment2, while sequence 65 has 138 gaps.
str_count(Cotesia.28S.alignment2[[1]], "-")
str_count(Cotesia.28S.alignment2[[65]], "-")

#Note that, when counted like this, gaps in flanking regions and internal gaps are both counted. Sometimes, we might rather trim gaps from the ends and then only count internal gaps or Ns.

#What if we want to look at the average number of gaps across all sequences in the alignment? Fortunately, again, our alignment (DNAStringSet object) is of a class that can be treated like a list. Therefore, we can use lapply().
lapply(Cotesia.28S.alignment1, str_count, ("-"))

#lapply() is an extremely helpful function that applies a function to every element of a list. Above, our first argument to lapply() is the alignment object (a DNAStringSet). The second argument is the function we want to apply to all elements. Note that the arguments for the function we want to apply (in this case str_count) need to be placed after a comma, within the arguments to be passed to lapply. So, please take note of this syntax if you are new to lapply. Check out the documentation for the base R function. Note that the line of code above is shorter than a loop would be to write. As well, "vectorized" functions also typically run faster. Conceptually, we can think of this basically as a loop in which we are doing the SAME ACTION on the elements of a list (or list-like object). This is a typical R way of doing things and is efficient to write and for computation.
?lapply

#If we want to look at a histogram showing the distribution of gap content across sequences. (See further commenting below about the unlist function.)
hist(unlist(lapply(Cotesia.28S.alignment1, str_count, ("-"))))

#The above is quite a nested way to write that. Let's tidyverse it up using piping! As well, here we are adding labels for our figure. Do you find this easier to read? Piping lets us organize our script better into a logical flow. Do this, then do that, then do that other thing. By contrast, with the nested coding structure above, we have to read the code from the inside out.
Cotesia.28S.alignment1 %>%
  lapply(str_count, ("-")) %>%
  unlist %>%
  hist(xlab = "Number of Gaps", ylab = "Frequency", main = "Distribution of Gap Count Among Cotesia 28S Sequences")

#What if we want to find the mean gap count? Note the following line won't work. You could uncomment it and try to run it as a demo. Using mean() as we've done so far won't work on our list object. Again, this is included as an example and will give an error. You can check it out by uncommenting and trying to run this line. Remember that we can check what type of data a function will accept by looking at its documentation. Check out the documentation for the base R function mean(). It is always wise to look up what kind of input data our function of choice is expecting.
#mean(lapply(Cotesia.28S.alignment1, str_count, ("-")))

#If we want to find the mean of all of the gap counts across list elements, we could do the following. The below is first using lapply() to apply the function str_count() to all of the elements in our alignment. We are counting the number of gaps for each sequence in the alignment. Then, we use unlist() to change the data type (from a list to a vector, a numeric vector in this case). This will enable us to pass the gap counts to the function mean() to calculate an overall mean. Note that the below doesn't save anything. Below, we are outputting the result to the screen. We could readily add a name and an assignment operator at the beginning if we wanted to save the results for further analysis.
mean(unlist(lapply(Cotesia.28S.alignment1, str_count, ("-"))))

#The above is another example of nested functions, written in base R style. In order to read that, remember that you need to start reading from the most nested set of parentheses and read outward. The order of running the functions will proceed from the inside out. If you have 2-3 functions, it is relatively straightfoward to read. However, here, we have four, as there is also a function as an argument to lapply. As well, you can imagine that this nesting could easily get out of hand if you want to add more steps. As well, it gets even worse to read if you need to specify multiple cases of non-default arguments for various functions. If you would like further reading on the topic of pipes, I suggest to check out the chapter on pipes in "R for Data Science": https://r4ds.had.co.nz/pipes.html

#So, again, let's tidyverse it up and compare against the base R style. The following is a way to do the same thing in tidyverse style, making use of piping. Many people would consider this easier to read, as the flow is more linear. Note that below I didn't assign this result to any new object, but that could readily be added if you want to save the result.
Cotesia.28S.alignment1 %>%
  lapply(str_count, ("-")) %>%
  unlist %>%
  mean

#Repeating the calculation of the mean number of gaps, for alignment 2
mean(unlist(lapply(Cotesia.28S.alignment2, str_count, ("-"))))

#Same thing, tidyverse style. Note that tidyverse style code tends to be longer (more lines of code) but narrower (fewer characters per line of code), and more linear.
Cotesia.28S.alignment2 %>%
  lapply(str_count, ("-")) %>%
  unlist %>%
  mean

#We could instead find the median value, as we aren't sure of the expected distribution of gap counts. The following is in base R style. For practice, you may convert the below to tidyverse style.
median(unlist(lapply(Cotesia.28S.alignment1, str_count, ("-"))))
median(unlist(lapply(Cotesia.28S.alignment2, str_count, ("-"))))

#We could also look at the distributions of gap numbers between these two alignments.
hist(unlist(lapply(Cotesia.28S.alignment1, str_count, ("-"))))
hist(unlist(lapply(Cotesia.28S.alignment2, str_count, ("-"))))

#These are simple plots using the hist() function. I suggest to review chapters 1 and 5 in "R for Data Science". For example, we could plot these on the same plot (using colour to differentiate the two sets of gap counts) or as side-by-side plots.

#So, we now have an idea about how changing the gap penalty influenced some features of the resulting alignments. Setting a heavier gap penalty resulted in a shorter, less-gapped alignment. We would also expect increasing nucleotide dissimilarity at certain positions in the alignment. We could examine this additional expectation through looking at the distribution of genetic distances among sequences under different alignment scenarios.

#Discussion and future work: But, how would we choose a suitable gap penalty? What criterion would we use to make our choice? How would we choose which alignment is better from an evolutionary point of view? Those are good questions, and there isn't necessarily an easy answer. For a ribosomal RNA gene region such as 28S, we could use alignment methods informed by secondary structure models. If you are working with ribosomal RNA genes (rRNA), I suggest to check out the alignment tools available through silva as well as the R package R4RNA, available through the Bioconductor repository. Alignment could be especially important for building hypotheses about phylogenetic relationships among distant relatives, as determining homology is more difficult in the case of evolutionarily distant lineages. For closely related groups and in conserved genes, standard alignment methods tend to work adequately as there are relatively few insertions and deletions (indels) among close relatives, and hence there is a higher certainty of homology and fewer gaps in the alignment. For protein-coding genes, we could use the amino acid sequences rather than nucleotides alone to build a more biologically-informed alignment. Amino acid sequences are more evolutionarily conserved than nucleotide sequences, and hence it is more feasible to find regions of homology among more distant relatives. The amino acid substitution matrices are informed by both the mutability and frequency of the 20 amino acids. Note that if our starting point is nucleotide data (DNA or RNA), we need to find the correct reading frame for translating our sequences.

#Tip for DNA barcode data specifically: In a previous study, our research group has found that a gap penalty of -3000 for muscle yields biologically realistic alignments for COI DNA barcodes (e.g. real indels are recovered in multiples of 3 nucleotides, indicating an animo acid insertion or deletion).

#We can also examine whether our final biological conclusions are robust to a range of reasonable choices of parameters settings. That would be termed a type of "sensitivity analysis". How sensitive are our findings to alignment parameter choice?

#FURTHER WORK: Check out other alignment algorithms available through the msa and DECIPHER packages. Compare the impact of alignment parameters settings on COI vs. 28S.

#FURTHER WORK: Export a visualization of the alignment. HINT: You could try the msaPrettyPrint() function in msa.

####PART 7 - CLUSTERING AND CONSEQUENCES OF ANALYSIS CHOICES ON NUMBER OF MOTU----

#Here, we will explore the consequences of alignment and clustering choices for downstream biological findings. We will continue to work with the 28S data. I suggest also to explore the COI data.

#Converting data type to DNAbin for further downstream analysis. Specifically, we are converting the object class so that we can create a distance matrix and then cluster the sequences.
dnaBin.Cotesia.28S.1 <- as.DNAbin(Cotesia.28S.alignment1)

?as.DNAbin

#The distance matrix being generated in the below example is using the TN93 model of sequence evolution. Model settings would require consideration for your own study. There are models of varying complexity that consider biological realism in calculating the distances (e.g. consider unequal nucleotide frequencies). Pairwise deletion of missing data and gapped sites would normally be preferred when missing data are spread throughout the alignment, when considering all sequences. By contrast, you might choose complete deletion if you had an alignment with little missing data but a region of high gap content and uncertainty regarding homology in a specific section of the alignment. Note that if missing data are widely spread throughout the alignment, the "complete deletion" option would result in a lot of data being thrown away. I don't recommend complete deletion when you have scattered missing data in various alignment positions; pairwise deletion is often more appropriate for large datasets of conserved marker genes.

distanceMatrix1 <- dist.dna(dnaBin.Cotesia.28S.1, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#documentation for the dist.dna function from the ape package
?dist.dna

#Note that for a phylogenetic analysis, we would typically want to perform modeltesting. This is a procedure for selecting which model of nucleotide evolution is most suitable for the data set. In some cases, we might NOT want to do model testing. For example, we might use a commonly-used metric of genetic distance if we want to compare genetic divergence results directly with other studies. The K2P model of sequence evolution is very commonly used in DNA barcoding studies. By contrast, the BIN algorithm is based upon p-distances. Therefore, a deliberate choice of model should be made, with a brief justification provided. You could check out this function from the phangorn package regarding modeltesting in R:
#https://www.rdocumentation.org/packages/phangorn/versions/2.4.0/topics/modelTest

#For now, we are moving on using the TN93 model. The TN93 model allows for two different rates for the two types of transitions (A <-> G and C <-> T) and a different rate for transversions (i.e. purine <-> pyrimidine).

#Next, we will perform clustering of our sequences according to the single-linkage method and using a 0.02 (2%) divergence threshold. You can play with the options and read the documentation. The below is using the IdClusters() function from the package DECIPHER. I like that function from the DECIPHER package as it offers a number of methods as potential arguments. Single linkage is a simple agglomerative clustering method (e.g. it starts with all points separate and then starts joining them). The threshold indicates that a point will join a group if it is less than that distance from ANY member of an existing cluster. This can result in chaining, where some members of a cluster are quite disparate. The BIN algorithm, currently only available through BOLD for COI sequences, uses a refined version of single linkage (Ratnasingham and Hebert 2013. PLoS ONE). As well, you can output a visualization in the same step if you wish (by setting type to both). The dendrogram is visualized with colours representing the clusters. Note that the colours will repeat for large dendrograms.
clusters.Cotesia.28S <- IdClusters(distanceMatrix1,
                                   method = "single",
                                   cutoff= 0.02,
                                   showPlot = TRUE,
                                   type = "both",
                                   verbose = TRUE)

#We can check out the documentation:
?IdClusters

#Viewing the result. Note that above we asked for both the clusters and the dendrogram by setting the type argument to "both". So, the output from the above is a list. The first element of the list is a dataframe, containing the cluster information. i.e. Which individual is assigned to which cluster? The second element in the list is the dendrogram. We could have asked for just the clusters or just the dendrogram if we had wanted to, but here we wanted both.
class(clusters.Cotesia.28S)
clusters.Cotesia.28S

#Viewing just the clusters (first element of the list). So, the cluster numbers could be used as Molecular Operational Taxonomic Units for biodiversity study. Note we would need double square brackets to dive into the element and retrieve the data for analysis.
clusters.Cotesia.28S[1]
clusters.Cotesia.28S[2]
clusters.Cotesia.28S[[2]]
str(clusters.Cotesia.28S[2])

#I suggest to read the R formal language definition regarding indexing. And, here is a small tutorial on lists: http://www.r-tutor.com/r-introduction/list
  
#Below, we want to look into the first element of our list (which is a data frame containing the cluster information). Within that, we want to look at the first column (which is the only column in this case).
  
#Counting the number of unique clusters that are assigned to the specimens. In the above case, we got 8 clusters. The double square brackets indicate that we retrieving the first element of the list (i.e. the dataframe with the cluster information). The single set of single square brackets then indicate we are looking at the first column within that dataframe. 
length(unique(unlist(clusters.Cotesia.28S[[1]][1])))

#Next, let's try a range of clustering thresholds. How do these impact the number of clusters (Molecular Operational Taxonomic Units) recovered from the analysis? For example, how many clusters do we get if we set the clustering threshold to a very low level, such as 0.005? What about a higher level such as 0.1? What about 0.2? What if we want to look across a large range?

#For a brief introduction to loops, you might check out:
#https://www.r-bloggers.com/how-to-write-the-first-for-loop-in-r/

#There was also a swirl tutorial on "for loops". You could also check out other online tutorials. For example, I found the following helpful: Further information on loops and practices to avoid. For example, it is recommended to avoid using cbind() with each iteration. It saves computation time to set the dimensions of your output from the beginning, rather than growing the size of the ouput object with every iteraction.
#https://swcarpentry.github.io/r-novice-inflammation/15-supp-loops-in-depth/
#That post clarifies practices to avoid, such as using cbind() with each iteration.

#Here is one potential solution to this iteration challenge. Is your solution similar or different? I would be interested to see how you solved this!

#First, we can will setup a vector with our desired cutoff values for clustering. In this case, we will try 0 (no different between sequences) to 0.5 (sequences differ by a distance of 0.5, which would equate to 50% of positions in the case of p-distances). Then printing to screen to check.
cut <- c(seq(0, 0.5, by = 0.005))
cut
class(cut)

#So, we have a numeric vector consisting of the clustering thresholds we want to try in this example.

#We will now create a vector to house our results, then we will check the class and print to screen to check. In our case, it is helpful to generate an object as output, as we want to plot that object and may want to populate it further.
number.clusters <- vector("numeric", length(cut))
class(number.clusters)
number.clusters

#Here, we will create a second vector, to use for an alernative loop format, for comparison.
number.clusters1 <- vector("numeric", length(cut))

#Next, we will iterate the clustering of our dataset over our cutoff values using a "for loop". Note that I made the first line below generic. We could wish to change the vector cut above in the future, such as to change the upper bound or the increment size for the stepping the clustering thresholds. Then, the below loop would still be fine because the number of iterations is tied to the length of the vector cut. We then use the Idclusters function to cluster the data, using each value in the vector cut as a clustering threshold. We then count the number of unique clusters generated in each iteration.

#Notes on syntax: The [, 1] means that we are looking at the first column of the dataframe generated by each iteration of the clustering. Note that for writing such a loop, it works best to work one step at a time. For example, we could first check what class of object the IdClusters function creates when run using the type  argument set to "clusters". It creates a data frame. We can readily work with that and don't need to change the data class. We then assign the number of clusters generated through each iteration to our vector called number.clusters, which we had already created. We fill the vector number.clusters in sequence using the number of clusters generated at each clustering threshold.

for (i in 1:length(cut)) {
  clusters.temp <- IdClusters(distanceMatrix1, method = "single", cutoff = cut[i], type = "clusters")
  number.clusters[i] <- length(unique(clusters.temp[, 1]))
}

#Have a look at the output
number.clusters

#We can plot the number of clusters generated for each clustering threshold value. Cut should be our independent variable.
plot(cut, number.clusters)

#In this figure, we can see that the sequences become more lumped into fewer MOTU as we increase the clustering threshold.

#The for loop is a frequently used tool. However, the while loop is also suitable here, and I found this to be faster in this case even though a condition is checked with each iteration. The while loop is considered more general than the for loop.

#initializing i at 1
i <- 1

#while loop. Up to the length of our vector named cut (which contains our candidate clustering thresholds), the loop will be repeated. At the end, i is stepped by 1. Otherwise, this loop is the same as the for loop above.

while(i <= length(cut)) {
  clusters.temp1 <- IdClusters(distanceMatrix1, method = "single", cutoff = cut[i], type = "clusters")
  number.clusters1[i] <- length(unique(clusters.temp1[, 1]))
  i = i + 1
}

plot(cut, number.clusters1)

#Checking if the results are the same between the for loop and while loop. YES.
all.equal(number.clusters, number.clusters1)

#Humorous cartoon about the while loop: https://www.reddit.com/r/ProgrammerHumor/comments/a5mghb/the_importance_of_knowing_how_to_correctly_use/

#Timing the two approaches. Here, the loops are passed as arguments to the function system.time(). TIP: For more extensive benchmarking options, I suggest to check out the packages rbenchmark and microbenchmark. This posting compares several benchmarking tools: https://www.alexejgossmann.com/benchmarking_r/

system.time(while(i <= length(cut)) {
  clusters.temp1 <- IdClusters(distanceMatrix1, method = "single", cutoff = cut[i], type = "clusters")
  number.clusters1[i] <- length(unique(clusters.temp1[, 1]))
  i = i + 1
})

system.time(for (i in 1:length(cut)) {
  clusters.temp <- IdClusters(distanceMatrix1, method = "single", cutoff = cut[i], type = "clusters")
  number.clusters[i] <- length(unique(clusters.temp[, 1]))
})

#so, in this case, the while loop is faster.

#Further work: What alternative solutions did you come up with for the iteration problem?

#Further work: Explore the impact of clustering thresholds on COI data in Cotesia, as contrasted wtih 28S. You might wish to sample down the COI dataset for exploratory analysis. HINT: The "Data Wrangling with dplyr and tidyr" cheat sheet shows helpful functions for data wrangling and filtering. https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf


####PART 8 - EXAMPLE OF BUILDING PHYLOGENETIC TREES IN R: NEIGHBOUR-JOINING APPLIED TO COI----

#Above, we built a dendrogram. We can also use R to build more types of trees. We will start from scratch and explore the COI data instead this time. 

#We will read in the original data and then filter for the COI-5P marker and for the presence of a sequence. We will also filter out records lacking a BIN. Sequences bearing a BIN (a Barcode Index Number) are at least 500 bp long and with fewer than 1% Ns. So, this is a quality control mechanism. By then sampling one sequence per BIN, we will also have a sample that spans the phylogenetic diversity of the data set.

#Reading in original data download from BOLD.
Cotesia <- read_tsv("Cotesia_BOLD_data.tsv")

#Filtering the data as describe above.
Cotesia.COI.subset <- Cotesia %>%
  filter(markercode == "COI-5P") %>%
  filter(str_detect(nucleotides, "[ACGT]")) %>%
  filter(!is.na(bin_uri)) %>%
  group_by(bin_uri) %>%
  sample_n(1)

#A commonly-used and fast technique for building a phenogram is neighbour joining. nj() is a function from the ape package that will perform a neighbour-joining analysis on the provided distance matrix.

#First, converting the sequences to a DNAStringSet
Cotesia.COI.subset$nucleotides <- DNAStringSet(Cotesia.COI.subset$nucleotides)

#Aligning the sequences using the muscle algorithm and also converting to the DNAbin data type. As this is a multiple sequence alignment, I used msa at the beginning of the object name.
msa.Cotesia.COI.subset <- muscle::muscle(Cotesia.COI.subset$nucleotides)

#Adding sequence names to the alignment. Using the processid, which is a unique identifier.
rownames(msa.Cotesia.COI.subset) <- Cotesia.COI.subset$processid

#Writing to a file.
writeXStringSet(DNAStringSet(msa.Cotesia.COI.subset), file = "msaCotesiaCOI.fas", format = "fasta")

#Calculating a distance matrix. Note again that the function dist.dna wants a DNAbin object, and so I embedded a conversion function.
distanceMatrix.COI <- dist.dna(as.DNAbin(msa.Cotesia.COI.subset), model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#Checking the distribution of genetic distances. Checking whether there are any apparent outliers that should be explored further. I find the shape of the distribution as well as the range of values (up to 0.15 or 15%) to be within the expected range for a genus at the COI marker.
hist(distanceMatrix.COI)

# Now passing the distance matrix to nj(), which is a function from the ape package that will perform a neighbour joining analysis. Neighbour joining is a fast algorithm and is preferable to UPGMA in terms of recovering relationships in the face of variability in rates of molecular evolution. The original neighbour joining paper is: Saitou and Nei 1987. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Mol Biol Evol. 4(4):406-25. (The neighbour joining algorithm is also one of the optins for IdClusters from DECIPHER.)
COItree <- nj(distanceMatrix.COI)

#plotting the NJ phenogram
plot(COItree)

#same as:
plot.phylo(COItree)

#That figure is hard to read. Let's look at the options:
?plot.phylo

#Let's try a visualization without tip labels and in a fan shape. It's a little easier to see the shape of the tree. 
plot.phylo(COItree, show.tip.label = FALSE, type = "fan")

#For further visualization ideas, I suggest to check out the following package and vignette:
#https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html

#Can you come up with a cool visualization? Examples might include:
#Add a dot to the tips, coloured by the dominant country for each BIN.
#Reduce the dataset to records having both 28S and COI. Build a tree for each and then compare them using the tanglegram function from the package dendextend.
#Please share!

#Also, check out the package phangorn for other options for building trees in R, including maximum parsimony and maximum likelihood. ML is generally expected to outperform MP due to evidence of variability in rates of molecular evolution in real data. MP can suffer from problems such as long branch attraction. https://cran.r-project.org/web/packages/phangorn/vignettes/Trees.pdf


