##input file created by featureCount
##convert the featureCount output file to matrix format
##matrix format like
##  gene_name   sample1 sample2 sample3 ...
##    gene1       1        2      1     ...
##    gene2       9        7      8     ...
##     ...        ...     ...    ...    ...




#!/usr/bin/perl
@infiles=@ARGV;
my %gene_matrix,@allsamples;
foreach $infile (@infiles){
    open(IN,$infile) or die "$!\n";
    <IN>;
    $line=<IN>;
    chomp $line;
    @strs=split(/\s+/,$line);
    @samples=@strs[6..$#strs];
    for($i=0;$i<=$#samples;$i++){
        $samples[$i] =~ s/[\/\\\.].+//;
    }
    push(@allsamples,@samples);
    while($line=<IN>){
        chomp $line;
        @strs=split(/\s+/,$line);
        if($strs[0] eq ""){
            next;
        }
        for($i=0;$i<=$#samples;$i++){
            $gene_matrix{$strs[0]}{$samples[$i]}=$strs[$i+6];
        }
    }
}
#@allsamples = sort(@allsamples)
print "gene\t".join("\t",@allsamples)."\n";
foreach $gene (sort keys(%gene_matrix)){
    print $gene;
    foreach $sample (@allsamples){
        if(exists($gene_matrix{$gene}{$sample})){
            print "\t".$gene_matrix{$gene}{$sample};
        }else{
            print "\t0";
        }
    }
    print "\n";
} 
