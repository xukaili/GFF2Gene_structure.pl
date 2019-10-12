#!/usr/bin/perl
my $usage=<<USAGE;
Usage:
     perl  $0  genome.gff  Your_gene_ID
e.g. perl  $0  test.gff    Si1g001100
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
        @gff = split(/\t/,$_);
        $mychr = $gff[0];
        if ($gff[6] eq '+') {
            $mystart = $gff[3] - 2000;
            $myend   = $gff[4] + 500;
        }
        elsif ($gff[6] eq '-') {
            $mystart = $gff[3] - 500;
            $myend   = $gff[4] + 2000;
        }
    }
}
close GFF;

open OUTR,">Plot_gene-$ARGV[1].R" or die $!;
print OUTR <<EOF;
genemodel_plot<-function(model, xaxis=TRUE) {
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
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-2010, end+2500), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="CDS") {
        rect(model\$newstart[i], 0, model\$newstop[i], .2, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="mRNA") {
        arrows(x0=model\$newstop[i]+2000,y0=.1,x1=model\$newstart[i]-500,y1=.1, lwd=1, col="dodgerblue4", length = 0.1)
        #arrows(x0=model\$newstop[i]+.2*length,y0=.1,x1=model\$newstart[i]-.2*length,y1=.1, lwd=1, col="dodgerblue4", length = 0.1)
      }
      if (type=="three_prime_UTR") {
        x<-c(model\$newstop[i], model\$newstop[i], model\$newstart[i]+.02*length, start, model\$newstart[i]+.02*length)
        y<-c(0,.2,.2,.1,0)
        polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0, model\$newstop[i], .2, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (orientation=="+") {
    model\$newstart<-start+model\$start-1
    model\$newstop<-start+model\$bpstop-1
    model<-model[which(model\$feature!="exon"),]
    model<-model[which(model\$feature!="gene"),]
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-2500, end+2010), ylim=c(-9, 0.4))
    for (i in 1:nrow(model)) {
      type<-model\$feature[i]
      if (type=="mRNA") {
        arrows(x0=model\$newstart[i]-2000,y0=.1,x1=model\$newstop[i]+500,y1=.1, lwd=1, col="dodgerblue4", length = 0.1)
        #segments(x0=model\$newstart[i],y0=.1,x1=model\$newstop[i],y1=.1, lwd=1, col="dodgerblue4")
      }
      if (type=="CDS") {
        rect(model\$newstart[i], 0, model\$newstop[i], .2, col = "steelblue3", border="dodgerblue4", lwd=1)
      }
      if (type=="three_prime_UTR") {
      	x<-c(model\$newstart[i], model\$newstart[i], model\$newstop[i]-.02*length, end, model\$newstop[i]-.02*length)
      	y<-c(0,.2,.2,.1,0)
      	polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=1)
      }
      if (type=="five_prime_UTR") {
        rect(model\$newstart[i], 0, model\$newstop[i], .2, col = "lightsteelblue1", border="dodgerblue4", lwd=1)
      }
    }
  }
  if (xaxis==T)   {
    Axis(side=3, labels=T, cex.axis=0.7)
  }
  rect(start -1500, -0.5, start -1000, -0.7, border = "dodgerblue4" , col ="lightsteelblue1")
  text(start -1600, -0.6, "UTR", cex=0.7, col="black", pos=4, offset=-1)
  rect(start -1500, -1, start -1000, -1.2, border = "dodgerblue4" , col ="steelblue3")
  text(start -1600, -1.1, "CDS", cex=0.7, col="black", pos=4, offset=-1)

}


pdf(file="Plot_gene-$ARGV[1].pdf")
$ARGV[1] <- read.table('gene_gff-$ARGV[1].txt',stringsAsFactors = F, header = F,comment.char = "#",sep = '\\t')
genemodel_plot(model=$ARGV[1], xaxis=T)

dev.off()
EOF


system("Rscript Plot_gene-$ARGV[1].R");
system("rm -rf gene_gff-$ARGV[1].txt");
system("rm -rf Plot_gene-$ARGV[1].R");

