require 'bio/appl/seala/blosum'

module Bio
  module Alignment
    class OriginalAlignment
      
      def karlin96
        blosum = Blosum.new open(File.join(File.dirname(`which calcd`),
                                           "blosum62.bla"))
        blosum.karlin!
        @score =
          valind.map do |ci|
            originalcolumn = self.slice(ci..ci).values
            column = originalcolumn.reject{|i| i=='.'}
            if column.size*2 < originalcolumn.size
              'nan'
            else
              n = column.size
              spscore = column[0..-2].zip((0...n-1).to_a).sum { |i|
                          column[i[1]+1..-1].sum { |j| blosum[i.first][j] }}
              (n>1) ? spscore / (n * (n-1) / 2).to_f : nil
            end
          end
      end
      
      def thompson97
        blosum = Blosum.new open(File.join(File.dirname(`which calcd`),
                                           "blosum62.bla"))
        aa = blosum.keys[0..-5].dup
        @score =
          valind.map do |ci|
            originalcolumn = self.slice(ci..ci).values
            column = originalcolumn.reject{|i| i=~/[.BXZ]/}
            if column.size*2 < originalcolumn.size
              'nan'
            else
              n = column.size
              vectors = column.map do |ci|
                          Vector[*aa.map{ |ai| blosum[ai][ci] }]
                        end
              centervector = vectors.sum.map{|vi| vi/n.to_f}
              #puts centervector; gets
              diffs = vectors.map{ |vi| (vi - centervector).r }
              #puts diffs; gets
              diffs.sum / n.to_f * (n.to_f / self.size)
            end
          end
      end
      
    end
  end
end
