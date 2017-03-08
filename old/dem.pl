#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

############################################################# PART I ########################################################################################
#Reads the parent directory ($parent) i.e./ifs/data/Hiseq/Runs recursively entering into each run folder, skipping if 'Unaligned' already exists,#
#checks presence of RTAcomplete.txt, checks presence of original samplesheet.csv, if yes then edits samplesheet to create new samplesheets as required.
# The hash, '%folders_to_be_processed' now contains runs with newly created sample sheets which needs to be demultiplexed.
#######################################################################################################################################################

my $par_dir; 
my %folders_to_be_processed;
my $parent = "/ifs/data/Hiseq/Runs";
opendir($par_dir,$parent);
my $path;
OUTER:
while (my $sub_dir = readdir($par_dir)){
       next if($sub_dir =~/^..?$/); 
       $path = $parent . '/' . $sub_dir;
       next unless(-d $path);
       opendir (my $sub_folder,$path);
       INNER:while (my $fol = readdir($sub_folder)){
                 next if($fol =~/^..?$/);
                 my $fol_path = $path . '/' . $fol;
                 my @Unaligned = "";
                 @Unaligned = glob("$path/Unaligned*");
                 next OUTER  if (@Unaligned);
                 my $RTAComplete = $path.'/'."RTAComplete.txt";
                 if (-f $RTAComplete){
                     if ($fol_path =~/\.csv$/){
                         $folders_to_be_processed{$path}++;
                         my $sample_sheet= $fol_path;
                        #system("dos2unix $sample_sheet");
                         system("sed -i 's/\r//g' $sample_sheet");
                         open (my $infile,"<",$sample_sheet)|| die ("cant open $!");
                         my %hash;
                         while (my $line = <$infile>) {
                                chomp $line;
                                next if ($line=~/^FCID/);
                                my @arr = split(/\,/,$line);
                                $arr[4]=~s/\s+//g;
                                $arr[5]=~s/\s+//g;
                                $arr[2]=~s/\s+/_/g;
                                $arr[2]=~s/\-/_/g;
                                $arr[2]=~s/\.//g;
                                $arr[2]=~s/\(//g;
                                $arr[2]=~s/\)//g;
                                if ($arr[5]=~/[ATGC]{4,}/){
                                    $arr[4]=$arr[4]."-".$arr[5];
                                    splice @arr,5,1;
                                    open (my $OUT1,">>","$path/SampleSheet_dualIndex.csv")|| die ("cannot open $!");
                                    print $OUT1 join(",",@arr),"\n";
                                    close($OUT1);
                                }    
                                else{
                                     my $index_length = length($arr[4]);
                                     $hash{$index_length}++;
                                     foreach my $index_length(keys %hash){
                                          if (length($arr[4])== $index_length){
                                              my $arr_len = @arr;
                                              if ($arr_len < 10){
                                                  my $diff =10 - $arr_len;
                                                  $arr[-1]= ("," x $diff).$arr[-1];
                                              }
                                              if ($arr_len >10){
                                                  my $diff =$arr_len -10;
                                                  splice @arr,5,$diff;
                                              }
                                              my $fh;
                                              open ($fh,">>","$path/SampleSheet_I$index_length.csv")|| die "cant open $!";
                                              print $fh join(",",@arr),"\n";
                                              close ($fh);
                                         }
                                     }
                                }
                         }#end of while block reading original samplesheet.csv file and creating new samplesheet.csv's
                     }#end of if block where its checks if samplesheet is added in run folder or not by Laurence 
                 }#end of if block after checking presence of RTAcomplete
             }#end of while loop that reads each run folder in the parent root directory
};#end of while loop reading root directory

############################################################# Part II #######################################################################
# For each run to be processed extract the information form 'RunInfo.xml' and accordingly create arguments for --use-bases option.
# Create job scripts with all commands and fire jobs simeltaneously.
#############################################################################################################################################

foreach my $process_folder(keys %folders_to_be_processed){
print $process_folder,"\n";
my $i1  ; my $i2 ; my $r2 ;
opendir (my $sub_folder,$process_folder);
while (my $files = readdir($sub_folder)){
       my $file = $process_folder."/".$files;
       if ($file=~/RunInfo.xml$/){
           my $xml = $file;
           ($i1,$i2,$r2) = parseXML($xml);
       }     
}close($sub_folder);
my $run = "";
my $index = "";
if (defined$i2){$index = "dual_index";} else {$index = "One_index";}
if (defined$r2){$run = "paired_end";} else {$run = "single_end";}
my @SSarray = glob("$process_folder/SampleSheet*");
OUTER2:
foreach my $csv (@SSarray){
        system("sed -i '1i FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,Project' $csv");
        open (my $csv_file,"<","$csv")|| die "cannot open $!";
        while (my $line = <$csv_file>){
               chomp $line;
               next if ($line=~/^FCID/);
               my @cols = split(/\,/,$line);
               my $n = "n"; 
               if ($cols[4]=~/[ATGC]-[ATGC]/){    ## If sample sheet is dual index
                   my $jobfile = $process_folder."/"."job_duInd.sh";
                   my ($a,$b) = split(/_/,basename($csv));
                   $b =~s/\.csv//;
                   my ($x,$y) = split(/-/,$cols[4]);
                   my $in_len = length($x);
                   if ($in_len < $i1){    
                       if (defined $r2){  
                            open (my $OUT,">","$jobfile")|| die "cannot open $!";
                            print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len$n\*,I$in_len$n\*,y\*n --output-dir $process_folder/Unaligned_$b\n\n";
                            print $OUT "make -C $process_folder/Unaligned_$b -j 8 all" ;
                            close($OUT);
                          #  system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");                        
                       }
                       else{
                            open (my $OUT,">","$jobfile")|| die "cannot open $!";
                            print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len$n\*,I$in_len$n\* --output-dir $process_folder/Unaligned_$b\n\n";
                            print $OUT "make -C $process_folder/Unaligned_$b -j 8 all";
                            close ($OUT);
                          # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                            system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                      }  
                   }    
                   if ($in_len == $i1){ 
                       if (defined $r2){ 
                            open (my $OUT,">","$jobfile")|| die "cannot open $!";
                            print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len,I$in_len,y\*n --output-dir $process_folder/Unaligned_$b\n\n";
                            print $OUT "make -C $process_folder/Unaligned_$b -j 8 all" ;
                            close($OUT);
                          # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                            system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                       }
                       else{
                            open (my $OUT,">","$jobfile")|| die "cannot open $!";
                            print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len,I$in_len --output-dir $process_folder/Unaligned_$b\n\n";
                            print $OUT "make -C $process_folder/Unaligned_$b -j 8 all";
                            close($OUT);
                         # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           }
                   }
               next OUTER2; # to make sure only one job script per samplesheet is created rather , skip to outer loop without reading further lines in samplesheet.
               }
               else{  #i.e col[4] does not have '-' character.                  
                   my ($a,$b) = split(/_/,basename($csv));
                   $b =~s/\.csv//;
                   my $jobfile = $process_folder."/"."job_$b.sh";
                   my $in_len = length($cols[4]);
                   if ($index eq "One_index"){
                       if ($in_len < $i1){    
                           if (defined $r2){  
                               open (my $OUT,">","$jobfile")|| die "cannot open $!";
                               print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len$n\*,y\*n --output-dir $process_folder/Unaligned_$b\n\n";
                               print $OUT "make -C $process_folder/Unaligned_$b -j 8 all" ;
                               close($OUT);
                           #   system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                            system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                             
                           }
                           else{
                                open (my $OUT,">","$jobfile")|| die "cannot open $!";
                                print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len$n\* --output-dir $process_folder/Unaligned_$b\n\n";
                                print $OUT "make -C $process_folder/Unaligned_$b -j 8 all";
                                close($OUT);
                           #   system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                               system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           }
                       }
                       if ($in_len == $i1){ 
                           if (defined $r2){ 
                                open (my $OUT,">","$jobfile")|| die "cannot open $!";
                                print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len,y\*n --output-dir $process_folder/Unaligned_$b\n\n";
                                print $OUT "make -C $process_folder/Unaligned_$b -j 8 all" ;
                                close($OUT);
                             #  system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                               system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           }
                           else{
                                open (my $OUT,">","$jobfile")|| die "cannot open $!";
                                print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len --output-dir $process_folder/Unaligned_$b\n\n";
                                print $OUT "make -C $process_folder/Unaligned_$b -j 8 all";
                                close($OUT);
                            # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                               system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                               }
                       }
                   }# "One_index" sequecing run loop complete
                   if ($index eq "dual_index"){
                       if ($in_len < $i1){    
                           if (defined $r2){  
                               open (my $OUT,">","$jobfile")|| die "cannot open $!";
                               print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len$n\*,n\*,y\*n --output-dir $process_folder/Unaligned_$b\n\n";
                               print $OUT "make -C $process_folder/Unaligned_$b -j 8 all" ;
                               close($OUT);
                             # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                              system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           }
                           else{
                                open (my $OUT,">","$jobfile")|| die "cannot open $!";
                                print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len$n\*,n\* --output-dir $process_folder/Unaligned_$b\n\n";
                                print $OUT "make -C $process_folder/Unaligned_$b -j 8 all";
                                close($OUT);
                              #  system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                                 system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           }
                       }
                       if ($in_len == $i1){ 
                           if (defined $r2){ 
                                open (my $OUT,">","$jobfile")|| die "cannot open $!";
                                print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len,n\*,y\*n --output-dir $process_folder/Unaligned_$b\n\n";
                                print $OUT "make -C $process_folder/Unaligned_$b -j 8 all" ;
                                close($OUT);
                               # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                                system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                           }
                           else{
                                open (my $OUT,">","$jobfile")|| die "cannot open $!";
                                print $OUT "/usr/local/bin/configureBclToFastq.pl --input-dir $process_folder/Data/Intensities/BaseCalls/ --sample-sheet $csv --fastq-cluster-count=0 --use-bases-mask y\*n,I$in_len,n\* --output-dir $process_folder/Unaligned_$b\n\n";
                                print $OUT "make -C $process_folder/Unaligned_$b -j 8 all";
                                close($OUT);
                              # system("qsub -cwd -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                                system("qsub -cwd -V -S /bin/bash -e $process_folder -o $process_folder -M harshal\.inamdar\@csc\.mrc\.ac\.uk -m bea $jobfile");
                               }
                       }
                       
                   }#dual_index sequencing run, but only one index in csv, loop complete
               next OUTER2;
               }# end loop :if col[4] in csv has no "-".
        }#end while loop ; i.e while opening each csv file
}#end foreach loop i.e for each sample sheet within process_folder
}#end foreach loop i.e for each process_folder

###############################################################################################################################################
##Subroutine to parse Runinfo.xml
#################################
sub parseXML {
    my $xml =$_[0];
    my ($i1,$i2,$r2);
    open (my $xmlFile,"<",$xml)|| die "cannot open $!";
    while (my $line =<$xmlFile>) {
       chomp $line;
       if ($line=~/Number\=\"2\"/ and $line =~/IsIndexedRead\=\"Y\"/ and $line =~/NumCycles\=\"(\d+)\"/){
              $i1 = $1;
           }
       if ($line=~/Number\=\"3\"/ and $line =~/IsIndexedRead\=\"Y\"/ and $line =~/NumCycles\=\"(\d+)\"/){
             $i2 = $1;
           }
       if ($line=~/Number\=\"3\"/ and $line =~/IsIndexedRead\=\"N\"/ and $line =~/NumCycles\=\"(\d+)\"/){
             $r2 = $1;
           }
       if ($line=~/Number\=\"4\"/ and $line =~/IsIndexedRead\=\"N\"/ and $line =~/NumCycles\=\"(\d+)\"/){
             $r2 = $1;
           }
    }
    return ($i1,$i2,$r2);
}
#############################################################################################################
exit;
