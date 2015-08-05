#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;

# Converts from maf to aln. Only takes blocks present in all sequences
# bzip2 -d -c sequences.maf.bz2 | ./maf2aln.pl --exclude samples -o
# sequences.aln
#

sub in_excluded($$)
{
   my ($sample, $excluded) = @_;

   my $skip = 0;
   foreach my $excluded_sample (@$excluded)
   {
      if ($excluded_sample eq $sample)
      {
         $skip = 1;
         last;
      }
   }

   return($skip);
}

#
# Main
#

my ($exclude, $output);
GetOptions ("exclude|e=s"  => \$exclude,
            "output|o=s"   => \$output)
   or die("Incorrect options provided\n");

my (@exclude_samples, @samples);
if (defined($exclude))
{
   @exclude_samples = split(",", $exclude);
}

my %sequences_out;

# Read from file, or STDIN
my $i = 1;
while (my $line_in = <>)
{
   chomp $line_in;

   # Process sample line
   if ($line_in =~ /^# hal/)
   {
      $line_in =~ m/^# hal \((.+)\)/;
      my @in_sequences = split(",", $1);

      # parse into sample names
      foreach my $in_sequence (@in_sequences)
      {
         $in_sequence =~ m/^(.+):\d+$/;
         my $name = $1;

         unless (in_excluded($name, \@exclude_samples))
         {
            push(@samples, $name);
         }
      }
      @samples = sort @samples;

   }
   # Process alignment blocks
   elsif ($line_in =~ /^a/)
   {
      if ($i % 1000 == 0)
      {
         print STDERR "$i alignments processed\n";
      }
      $i++;

      my %found_seqs;
      my $seq_length;

      while (my $seq_line = <>)
      {
         chomp $seq_line;

         # Sequence in alignment block
         if ($seq_line =~ /^s/)
         {
            my ($s, $name, $num1, $num2, $strand, $pos, $seq) = split("\t", $seq_line);

            $name =~ m/^(.+?)\./;
            $seq_length = length($seq);

            unless(in_excluded($1, \@exclude_samples) || $1 eq "Anc0")
            {
               $found_seqs{$1} = $seq;
            }
         }
         # End of alignment block, process previous s lines
         else
         {
            # Check min number of samples are present
            if (scalar(keys %found_seqs) >= 0.2*scalar(@samples))
            {
               foreach my $sample (@samples)
               {
                  if (defined($found_seqs{$sample}))
                  {
                     $sequences_out{$sample} .= $found_seqs{$sample};
                  }
                  else
                  {
                     $sequences_out{$sample} .= "-" x $seq_length;
                  }
               }
            }

            last;
         }
      }
   }
}

# Open output fastas, and write sequence to them
my @tmp_files;
foreach my $sample (@samples)
{
   my $tmp_file = "tmp_$sample.fa";
   push(@tmp_files, $tmp_file);

   my $file_out = Bio::SeqIO->new( -file   => ">$tmp_file",
                                   -format => "fasta" ) || die ($!);

   $file_out->write_seq(Bio::Seq->new(-display_id => $sample, -seq =>$sequences_out{$sample}));
}

# cat files together, delete tmps
my $cat_command = "cat " . join(" ", @tmp_files) . " > $output";
system($cat_command);

my $rm_command = "rm " . join(" ", @tmp_files);
system($rm_command);

exit(0);

