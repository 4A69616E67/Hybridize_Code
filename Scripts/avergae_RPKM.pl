##usage: perl averagr_RPKM.pl <bedgraph file> <chromosome region>
##example: perl averagr_RPKM.pl TeLs-s10.bedgraph Chr3:153,637,230-153,637,330

use List::Util qw(min max);
$in_file = $ARGV[0];
$region=$ARGV[1] ;
$region =~ s/,//g ;
@strs = split(/:|-/,$region);
$chr=$strs[0];
$start=$strs[1];
$end=$strs[2];
open(IN,$in_file);
$flag=0;
$sum=0;
while($line=<IN>){
	chomp $line;
	@strs=split(/\s+/,$line);
	if($strs[0] eq $chr){
		if($strs[1]<$end && $strs[2]>$start){
			# print "$line\n";
			$sum+=$strs[3]*(min($end,$strs[2])-max($start,$strs[1]));
			$flag=1;
			next;
		}
	}
	if($flag==1){
		last;
	}
}
$res=$sum/($end-$start);
print "$chr\t$start\t$end\t$res\n";