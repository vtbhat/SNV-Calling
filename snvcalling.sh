#!/bin/bash

realign=0 #Save status of realignment
gunzip=0 #Save status of zipping VCF file
v=0 #Save status of verbose mode
index=0 #Save status of indexing
answer=0 #Check whether the input file need to be overwritten
while getopts "a:b:r:eo:f:zvijh" opt
do
	case $opt in
		a) reads1=$OPTARG;;
		b) reads2=$OPTARG;;
		r) ref=$OPTARG;;
		e) realign=1;;
		o) output=$OPTARG;;
		f) millsFile=$OPTARG;;
		z) gunzip=1;;
		v) v=1;;
		i) index=1;;
		j) answer=1;;
		h) echo "-a Location(reads1) -b Location(reads2) -r reference -o outputfilename -e -f millsFile -z -v -i -j"
	esac
done

if [ $v -eq 1 ]
then
	echo "Checking if all the files exist"
fi

#Check if reads1 exists
if [ ! -f "$reads1" ]
 then
 	echo "$reads1 does not exist"
 	exit
fi

#Check if reads2 exists
if [ ! -f "$reads2" ]
 then
 	echo "$reads2 does not exist"
 	exit
fi

#Check if ref exists
if [ ! -f "$ref" ]
 then
 	echo "$ref does not exist"
 	exit
fi

#Check if the output file exists, or must be overwritten
if [ -f "$output.vcf" ] || [ -f "$output.vcf.gz" ]
 then
 	echo "The output file exists. Do you wish to overwrite it? (Y or n)"
 	read reply
 	case reply in
 		Yes|Y|y) continue;;
 		No|N|n) exit;;
 	esac
fi

#Check if reads_1 are zipped
zipfilename=$(basename $reads1) 
zipext=${zipfilename##*.}
if [ $zipext == "gz" ]	
then
 gunzip -d -k $reads1
 reads1=$(sed 's/.gz//' <<<$reads1)
fi

#Check if reads are zipped
zipfilename=$(basename $reads2) 
zipext=${zipfilename##*.}
if [ $zipext == "gz" ]	
then
 gunzip -d -k $reads2 
 reads2=$(sed 's/.gz//' <<<$reads2)
fi


#Index the reference
if [ $v -eq 1 ]
then 
	echo "Indexing the reference"
fi
bwa index $ref




#Map reads to the reference
if [ $v -eq 1 ]
then 
	echo "Mapping the reads to the reference"
fi
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $ref $reads1 $reads2 > aligned.sam

#If answer=ON, replace FASTQ with SAM
if [ $answer -eq 1 ] && [ -f "aligned.sam" ]
then
	rm $reads1 
	rm $reads2
fi


#Convert to bam and sort the aligned file
if [ $v -eq 1 ]
then 
	echo "Converting to BAM and sorting the aligned file"
fi

samtools sort -O bam -o aligned.bam -T /tmp aligned.sam

#If answer=ON, replace FASTQ with SAM
if [ $answer -eq 1 ] && [ -f "aligned.bam" ]
then
	rm "aligned.sam"
fi


##Improvement of the alignment
#For the realignment of BAM file
if [ $realign -eq 1 ]
then
	if [ $v -eq 1 ]
	then 
	echo "Improving the alignment: creating .faidx and .dict files"
	fi
	#Index reference sequece
	samtools faidx $ref
	#Find the basename of the FASTA reference file
	dictfilename=$(basename $ref .fa) 
	dictname="$dictfilename".dict
	#Create a sequence dictionary
	samtools dict $ref -o $dictname
	samtools index aligned.bam
if [ $v -eq 1 ]
	then 
	echo "Improving the alignment: performing realignment. Ensure JAVA 8 is installed and the default java (required for GATK 3.7) and the gatk folder is in your PATH"
	fi
	java -Xmx2g -jar gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I aligned.bam -o aligned.intervals -known $millsFile &> $output.log
	java -Xmx4g -jar gatk/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I aligned.bam -targetIntervals aligned.intervals -known $millsFile -o realigned.bam &>> $output.log
fi

#For indexing of BAM file
if [ $index -eq 1 ]
then
	if [ $v -eq 1 ]
	then 
	echo "Indexing the BAM file"
	fi
	samtools index realigned.bam
fi


#Call variants
if [ $v -eq 1 ]
then 
	echo "Creating mpileup and calling the variants"
fi
bcftools mpileup -Ou -f $ref realigned.bam | bcftools call -vmO z -o $output.vcf


#Create BED file from VCF
if [ $v -eq 1 ]
then 
	echo "Creating the BED file from the VCF file $output"
fi
#Delete the header lines
sed '/^#/d' $output.vcf > temp.txt
#Obtain CHR, POS, REF and ALT columns and remove CHr
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"(length($5)-length($4))}' temp.txt | awk '{print $1"\t"$2"\t"($2+$5)"\t"$5}' | sed 's/chr//' > temp1.txt
#Place SNPs and indels in different files
awk '{if($4==0) print $0}' temp1.txt > snps.txt
awk '{if($4!=0) print $0}' temp1.txt > indels.txt
rm temp.txt
rm temp1.txt


#To zip the vcf output
if [ $gunzip -eq 1 ]
then 
	if [ $v -eq 1 ]
then 
	echo "Zipping the VCF file"
fi	
	gzip $output.vcf
fi

if [ $v -eq 1 ]
then 
	echo "Completed the variant calling"
fi
#End of file
