#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  genome.gff  Your_gene_ID
e.g. perl   GFF2Gene_structure.pl   test.gff   Si1g08100
USAGE
print $usage if(@ARGV==0);
exit if(@ARGV==0);

$mychr = $mystart = $myend =0;
system("grep  '$ARGV[1]' $ARGV[0] > gene_gff-$ARGV[1].txt");

open GFF ,"gene_gff-$ARGV[1].txt" or die "$!";
while (<GFF>) {
    $_ =~ s/\n*|\r*//g;
    if($_ =~ /\texon\t/){
        @F = split(/\t/,$_);
        $F[8]=~/Parent=(.*?);/;
        $cds_pos{$1}{"$F[3] $F[4] $F[7]"}=1;
        $contig{$1}=$F[0];
        $stream{$1}=$F[6];
    }
    if ($_ =~ /\tgene\t/){
        $_ =~ /;Note=(\d+) /;
        $gene_n = $1;
        @gff = split(/\t/,$_);
        $mychr = $gff[0];
        if ($gff[6] eq '+') {
            $mystart = $gff[3] - 2000;
            $myend   = $gff[4] + 1000;
        }
        elsif ($gff[6] eq '-') {
            $mystart = $gff[3] - 1000;
            $myend   = $gff[4] + 2000;
        }
    }
}
close GFF;
$gene_n = `grep 'mRNA' gene_gff-$ARGV[1].txt |wc -l`;
$gene_n =~ s/\s+//g;
print "Query:  $ARGV[1]  $gene_n  :  $mystart  -  $myend\n\n";

for ($i =1; $i <= $gene_n; $i++) {
  $new_gffF = $ARGV[1] . '.' . $i . ';';
  $new_gff = $ARGV[1] . '.' . $i;
  system("grep  '$new_gffF' gene_gff-$ARGV[1].txt > tmp-gene_gff-$new_gff.txt");
}

open OUTR,">Plot_gene-$ARGV[1].R" or die $!;
print OUTR <<EOF;
genemodel_plot<-function(model, xaxis=TRUE, drop=0) {
  par(mar=c(1,1,3,1), cex=1)
  colnames(model) =c('chr','db','feature','start','end','tmp1','orientation','tmp3','gene')
  orientation = model\$orientation[1]
  start <- min(c(model\$start,model\$end))
  end <- max(c(model\$start,model\$end))
  tmp_min=min(c(model\$start,model\$end))
  model\$start=model\$start-tmp_min+1
  model\$end=model\$end-tmp_min+1
  tmp_max=max(c(model\$start,model\$end))
  tmp_min=min(c(model\$start,model\$end))
  model<-cbind(as.character(model[,3]), as.numeric(model[,4]), as.numeric(model[,5]) )
  model<-as.data.frame(model)
  colnames(model)<-c("feature", "start", "bpstop")
  model\$start<-as.numeric(as.character(model\$start));model\$bpstop<-as.numeric(as.character(model\$bpstop))
  length<-tmp_max-tmp_min
  if (orientation=="-") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c($mystart, $myend), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="CDS") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="mRNA") {
        arrows(x0=model\$newstop[i]+2000,y0=0.1 - drop,x1=model\$newstart[i]-1000,y1=0.1 - drop, lwd=1, col="dodgerblue4", length = 0.1)
      }
      if (type=="three_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (orientation=="+") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c($mystart, $myend), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="mRNA") {
        arrows(x0=model\$newstart[i]-2000,y0=0.1 - drop,x1=model\$newstop[i]+1000,y1=0.1 - drop, lwd=1, col="dodgerblue4", length = 0.1)
      }
      if (type=="CDS") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="three_prime_UTR") {
         rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0 - drop, model\$newstop[i], 0.2 - drop, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (xaxis==T)   {
    Axis(side=3, labels=T, cex.axis=0.7)
  }

}

pdf(file="Plot_gene-$ARGV[1].pdf")
EOF
for ($i = 1; $i <= $gene_n; $i++) {
  $new_gff = $ARGV[1] . '.' . $i;
  $j = ($i)*0.5 - 0.5;
  print OUTR "$new_gff <- read.table('tmp-gene_gff-$new_gff.txt',stringsAsFactors = F, header = F,comment.char = \"#\",sep = '\\t')\n";
  print OUTR "genemodel_plot(model=$new_gff, xaxis=T, drop=$j)\n";
  print OUTR "text($mystart + 1000, 0.3 - $j, \"$new_gff\", cex=0.7, col=\"black\", pos=4, offset=-1)\n";
  print OUTR "par(new = TRUE)\n";
}
  print OUTR "rect($mystart + 300, -0.5- $j, $mystart +800, -0.7- $j, border = \"dodgerblue4\" , col =\"lightsteelblue1\")\n";
  print OUTR "text($mystart + 200, -0.6- $j, \"UTR\", cex=0.7, col=\"black\", pos=4, offset=-1)\n";
  
  print OUTR "rect($mystart + 1500, -0.5- $j, $mystart +2000, -0.7- $j, border = \"dodgerblue4\" , col =\"steelblue3\")\n";
  print OUTR "text($mystart + 1400, -0.6 - $j, \"CDS\", cex=0.7, col=\"black\", pos=4, offset=-1)\n";
  
  print OUTR "text($mystart + 2500, -0.6- $j, \"$ARGV[1]\", cex=0.7, col=\"black\", pos=4, offset=-1)\n";
print OUTR "dev.off()\n";


system("Rscript Plot_gene-$ARGV[1].R");
system("rm -rf gene_gff-$ARGV[1].txt");
system("rm -rf Plot_gene-$ARGV[1].R");
system("rm -rf tmp*");

