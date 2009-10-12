require 'bio/appl/seala/blosum'

module Bio
  module Alignment
    class OriginalAlignment
      
      Liu08matrix =
'A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
10 4  3  3  6  4  4  6  3  4  4  4  4  3  4  7  6  2  3  6
4  10 5  3  2  6  5  3  5  2  3  7  4  2  3  4  4  2  3  2
3  5  10 6  3  5  5  5  6  3  3  5  3  3  3  6  5  2  3  3
3  3  6  10 3  5  7  4  4  3  2  4  3  3  4  5  4  2  3  3
4  3  3  3  10 3  2  3  3  4  4  3  4  3  3  4  4  3  3  4
4  6  5  5  2  10 7  3  5  2  3  6  5  2  4  5  4  3  4  3
5  5  5  7  2  7  10 4  5  3  3  6  4  3  5  5  5  3  4  4
5  3  5  4  3  3  3  10 3  2  2  3  3  3  3  5  3  3  3  3
2  4  5  3  2  4  4  2  10 2  2  3  2  3  2  3  2  2  5  2
5  3  3  3  5  3  3  2  3  10 8  3  7  6  3  4  5  3  5  9
5  4  3  2  5  4  3  2  3  8  10 4  8  6  3  4  5  4  5  7
4  7  5  4  2  6  6  3  4  2  3  10 4  2  4  5  4  2  3  3
4  4  3  2  4  5  3  2  3  6  7  4  10 5  3  4  4  4  4  6
3  3  3  3  3  3  3  3  4  5  5  3  5  10 2  3  3  6  8  4
4  3  3  4  2  4  4  3  3  2  2  4  3  2  10 4  4  2  2  3
7  4  7  6  4  6  6  6  4  3  3  6  4  3  4  10 7  2  3  3
4  3  4  3  3  3  3  2  2  3  3  3  3  2  3  6  10 2  2  4
2  2  2  2  3  3  2  3  3  2  3  2  3  4  2  2  3  10 5  2
3  3  3  2  3  3  3  2  6  3  3  3  3  7  2  3  3  6  10 3
6  2  2  2  4  3  3  2  2  9  7  3  7  4  3  3  6  2  4  10
'
      def liu08(w=1)
        blosum = Blosum.new Liu08matrix
        setweights! w
        valind.map do |ci|
          originalcolumn = self.slice(ci..ci).values
          column = originalcolumn.reject{ |i| !Alphabet.member? i}
          if column.size*2 < originalcolumn.size
            'nan'
          else
            n = column.size
            commonAA = column.mostcommon
            Alphabet.sum{|i|
              self.profile[ci][i] * blosum[commonAA][i] }
          end
        end
      end

      def liu08w
        liu08 3
      end

    end
  end
end

class Array
  def mostcommon
    keys = self.uniq
    histogram = keys.map{|k| self.count k}
    index = (0...histogram.size).zip(histogram).
              max{|h1,h2| h1[1]<=>h2[1]}.first
    keys[index]
  end
end
