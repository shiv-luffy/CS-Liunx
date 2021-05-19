## Linux command line basics    

### get current working directory    
```pwd```  

### list directory contents  
```ls```  
options:  
-l #gives a list with some detials  
-a #all file with hidden files  
-lh #gives list with human readable file sizes  

### Change/make/remove directory  
```cd```  
```mkdir dirName```  # creates a directory called dirName  
```rmdir dirName```  # removes a directory called dirName  
```.```  # present working directorys  
```..```  # directory above  
```~```  # user’s home directory  
```/```  # root directory  
```cd ../../home/test/```  #relative path  
```cd /home/test/```  #absolute path  
```|``` # pipe can used to combine commands
```>``` # redirect output into file
```>>``` # append output into file

### copy/move/remove  
```cp /path/to/fileName /path/to/destination/```  # copy thee file called fileName to the destination  
```mv /path/to/fileName /path/to/destination/```  # moves thee file called fileName to the destination  
```rm /path/to/fileName```  # removes the file called fileName  

### print something  
```echo "this is a test"``` #prints data  
```printf "this is a test \n"``` # formats and prints the data  

### Creating/editing files    
```nano fileName```  # creates if not exist or opens in editing mode the file called fileName  
```vim fileName```  # creates if not exist or opens in editing mode the file called fileName  
```vi fileName```  # creates if not exist or opens in editing mode the file called fileName  
```touch fileName```  # creates empty file called fileName  

###  Viewing a file  
```cat fileName``` # prints the contents of the file fileName on the terminal  
```zcat fileName.txt.gz # prints the contents of compressed file fileName on the terminal  
```less fileName```  # prints the contents of the file on window and not on the terminal  
```head fileName``` #prints the top 10 lines  
```tail fileName``` #prints the bottom 10 lines  
options:  
```-n <INT>``` #head and tail command takes this option and prints the ```<INT>``` number of line  

### counting
wc fileName``` # print lines, words, and byte counts for the file fileName  
options:  
```-w``` # print the word counts  
```-l``` # print the line counts  

### pattern matching/substitute
```grep 'pattern' fileName``` #grep searches the file fileName for lines containing a match to the given pattern  
options:  
```-o``` # print only the matched (non-empty) parts of a matching lines  
```-v``` # Invert matching  
```-c``` # print a count of matching lines  

```sed 's/oldPattern/newPattern/' fileName.txt``` # replaces the pattern “oldPattern” with “newPattern” in the file fileName  
```sed 's/oldPattern/newPattern/2' fileName.txt``` # replaces the pattern “oldPattern” with “newPattern” in the file fileName  
```sed 's/oldPattern/newPattern/g' fileName.txt``` # replaces the first, second occurrence of pattern “oldPattern” with “newPattern” in the file fileName  

### download files from Web  
```wget "www.domainName.com/test.txt"```  #downloads the text file  
```curl "www.domainName.com/test.txt"```  #downloads the text file  

### Compress/de-compress   
```tar -zcvf directory.tar.gz directory/```  # compress the directory  
```tar -zxvf archive.tar.gz```  # de-compress the archive  

### kill/suspend process, getting help  
```ctrl+c```  # kill a process  
```ctrl+z```  # suspending a process   
```man cmd```  # prints man page for the command ```cmd```  
```cmd --help```  # prints help on the teminal for the command ```cmd```  
```clear```  # clears the termianl window  


### Examples 

### Example 1  

##### small bash script   

```nano hello.sh```  #create bash file called hello  

##### enter these lines in the file
```
echo Hello $1  
echo How are you?  
```

##### now run the the scrpit with a user input
```bash hello.sh luffy```  

### Example 2   

##### print all the header line of fasta file  

```nano test.fasta```  #create text file called test.fasta  
##### enter these lines in the file
```
>test seq1  
atgacgatgctagctatctgactatgctgatc  
>test seq2  
atgacgatgctagctatctgactatgctgatc  
>test seq3  
atgacgatgctagctatctgactatgctgatc  
```
##### run these commands 
```grep ‘>’ test.fasta```  #prints all the headers (prints all lines with a pattern '>')  
```grep -c ‘>’ test.fasta```  #prints counts of lines with a pattern '>'  

### Example 3  

##### print all the header line of fasta file  

```nano test.vcf```  #create text file called test.vcf 
##### enter these lines in the file
```  
#This is a test file  
#some header info  
#some header info  
#CHROM	START	END	ID	REF	ALT  
chr1	18651	368746	.	G	C  
chr5	684338	8468468	.	C	T  
```  
##### run these commands   
```grep -v ‘#’ test.vcffffffffffff```  #prints all the lines not matching with ```#```  

