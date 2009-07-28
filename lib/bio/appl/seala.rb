require 'mathn'
require 'matrix'
require 'bio/appl/seala/lockless'
require 'bio/appl/seala/karlin'
require 'bio/appl/seala/entropy'
require 'bio/appl/seala/taylor'
require 'bio/appl/seala/liu'
require 'bio/appl/seala/cuttree'
require 'bio/appl/seala/fredrik-utilities.rb'

module Bio

  class Sequence
    class AA < String

      def aalength
        @aalength ||= self.gsub(/[^#{Bio::Alignment::Alphabet}]/,'').size
      end

      def ident(other)
        a = Bio::Alignment::Alphabet
        s = self.split('').zip(other.split('')).sum { |i|
              if i[0] && i[1] && (i[0]+i[1])=~/[#{a}]{2}/ && i[0]==i[1]
                1
              else
                0
              end }
        s.to_f / [self.aalength, other.aalength].max
      end

    end
  end

  module Alignment

    Alphabet = %w(A R N D C Q E G H I L K M F P S T W Y V)

    class OriginalAlignment

      attr_accessor :valind

      def validIndices
        return @valind if @valind
        @valind = []
          self[keys.first].split('').each_with_index do |ci,i|
            if ci !~ /\.|B|X|Z/
              @valind << i
            end
          end
        @valind
      end

      def distmatrix
        return @distmatrix if @distmatrix
        tmp = 
          self.keys.map do |seq1|
            self.keys.map do |seq2|
              (1 - self[seq1].ident(self[seq2]))
            end
          end
        @distmatrix = Matrix.rows tmp
      end

      def phylogeny
        return @phylogeny if @phylogeny
        L.mean [1] #Just a hack to load R
        @phylogeny =
          if self.size<2
            nil
          elsif self.size==2
            Rtree.new [[-1,-2]]
          else
            rmatrix = self.distmatrix.to_R
            Rtree.new R.agnes(R.as_dist(rmatrix),:diss=>'TRUE')['merge']
          end
      end

      def evoltrace(method = 'shannon')
        cache = Hash.new
        result =
          (1...self.size).sum do |node|
            g = self.phylogeny.groups(node)
            gsum = g.sum do |group|
                     group2 = if group.is_a? Fixnum
                                group-1
                              else
                                group.map{|gi|gi-1}
                              end
                     cache[group2] ||= begin
                                         al = self.subalignment(group2)
                                         al.valind = self.validIndices
                                         eval "al.#{method}"
                                       end
                   end.map{|i| i.to_f / g.size}
          end
        result.map{|i| i + 1}
        #result.map{|i| i.to_f / (self.size-1)}
      end

      # arr is a numeric array pointing to key indices
      def subalignment(arr)
        arr = [arr] if arr.is_a? Fixnum
        al = Bio::Alignment::OriginalAlignment.new
        arr.each do |i|
          next if (i>=self.size or i<0)
          key = self.keys[i].dup
          al.store(key,self[key].dup)
        end
        al
      end

      def alignfile
        if !@alignfile
          @alignfile = Tempfile.new('align')
          @alignfile.puts self.output_clustal
          @alignfile.flush
        end
        @alignfile
      end

      def weightfile
        if !defined? @weightfile
          tmp = self.w
          updateweightfile
        end
        @weightfile
      end

      def freqfile
        if !@freqfile
          setprofile!
        end
        @freqfile
      end

      def updateweightfile
        return if !defined? @w
        @weightfile = Tempfile.new('weight')
        (0...self.keys.size).each do |i|
          @weightfile.puts [i+1,keys[i],@w[keys[i]]].join("\t")
        end
        @weightfile.flush
        @weightfile
      end

      def w
        setweights! if !@w
        @w
      end

      def profile
        setprofile! if !@profile
        @profile
      end

      def setweights!(wmethod=1)
        @w = {}
        @profile = nil
        str = `calcw -i #{alignfile.path} -m #{wmethod}`
        str.split("\n").each do |line|
          @w[line.split[1]] = line.split.last.to_f
        end
        updateweightfile
        setprofile!
        @w
      end

      def setprofile!(assocmatrix=nil,beta=1,epsilon=0)
        param1 = Tempfile.new 'param1'
        param1.puts "50\n0 #{epsilon.to_f} #{beta.to_f}\n"
        param1.flush
        exec = "calcf -i #{alignfile.path} " + 
               "-w #{weightfile.path} -p #{param1.path} -t 1 "
        if assocmatrix
          raise "Can't find matrix file" if !File.exists? assocmatrix
          exec += "-c #{assocmatrix} -e "
        end
        output = `#{exec}`
        #puts '---',output,'---'
        @freqfile = Tempfile.new 'freq'
        @freqfile.puts output
        @freqfile.flush
        aa,*lines = output.split("\n")
        aa = aa.split[1..-1]
        @profile = []
        lines.each do |line|
          col,aai,p = line.split
          col = col.to_i
          aai = aai.to_i
          p   = p.to_f
          @profile[col] ||= Hash.new(0)
          @profile[col][aa[aai]] = p if p!=0
        end
        @profile
      end

      # Zvelebil counts properties
      def zvelebil87
        s = 
"Ile  Y             Y    
Leu  Y             Y    
Val  Y         Y   Y    
Cys  Y         Y        
Ala  Y         Y Y      
Gly  Y         Y Y      
Met  Y                  
Phe  Y               Y  
Tyr  Y     Y         Y  
Trp  Y     Y         Y  
His  Y Y   Y Y       Y  
Lys  Y Y   Y Y          
Arg    Y   Y Y          
Glu      Y Y Y          
Gln        Y            
Asp      Y Y Y Y        
Asn        Y   Y        
Ser        Y   Y Y      
Thr  Y     Y   Y        
Pro            Y       Y
#Asx        Y            
#Glx        Y            
#Gap  Y Y Y Y Y Y Y Y Y Y
#Unk  Y Y Y Y Y Y Y Y Y Y"
        #s = open('zvelebil').read.nsplit.reject{|i|i=~/#/}.inject({}) { |h,line|
        s = s.nsplit.reject{|i|i=~/#/}.inject({}) { |h,line|
              h[line.split.first.oneLetterAA] = line[5..-1]; h }
        self.validIndices.map do |ci|
          mem = self.slice(ci..ci).values.uniq.reject{ |i| !Alphabet.member? i}
          if (n=mem.size) == 1
            1
          else
            score = 0.9
            props = mem.map{|i| s[i].indices 'Y'}
            ind = props.flatten.uniq
            ind.each do |i|
              score -= 0.1 unless (props.count{|pi| pi.member? i} == n)
            end
            score
          end
        end
      end

      # Wu&Kabat is the simplest of scores.
      def wu70
        seala 1
      end

      def wu70w
        self.setweights! 3
        ans = seala 1
        self.setweights! 1
        ans
      end

      # Pei made the AL2CO software which has several methods though...
      def pei01var
        seala 2
      end

      def pei01sp
        seala 9
      end

      # Pei made the AL2CO software which has several methods though...
      def pei01varw
        self.setweights! 3
        ans = seala 2
        self.setweights! 1
        ans
      end

      def pei01spw
        self.setweights! 3
        ans = seala 9
        self.setweights! 1
        ans
      end

      # Sander and Schneider introduced the Shannon entropy
      # They use natural logarithm. Also Sum-of-Pairs.
      def sander91sp
        seala 10
      end
      def sander91shannon
        seala 3
      end

      # Shenkin introduced Shannon entropy at the same time as
      # Sander and Schneider. Shenkin used log2 and defined the final
      # score as 6 * 2^S.
      def shenkin91
        seala 4
      end

      # Gerstein and Altman measures entropy relative to entropy of a 
      # random alignment (entropy with a linear translation).
      def gerstein95
        seala 5
      end

      # Mirny uses Shannon entropy but with a smaller alphabet -
      # only six different residue types.
      def mirny99
        seala 6
      end

      # Williamson introduced Relative Entropy which is relative
      # to the amino acid distribution in the MSA. He also uses
      # a smaller alphabet of nine residue types.
      def williamson95
        seala 7
      end

      # Valdar and Thornton used a Sum-of-Pairs score.
      def valdar01
        self.setweights! 2
        ans = seala 8
        self.setweights! 1
        ans
      end

      # Valdar, in his review article, introduces a tricky score.
      # It is a product of entropy, stereochemical properties
      # and gaps. The problem is that the score is HEAVILY depending
      # on the weights for each of these factors.
      def valdar02
        seala 12
      end

      # Caffrey introduced the vonNeumann entropy.
      def caffrey04
        seala(11,'qij')
      end

      # Caffrey introduced the vonNeumann entropy.
      def caffrey04w
        self.setweights! 3
        ans = seala(11,'qij')
        self.setweights! 1
        ans
      end

      # rate4site
      def mayrose04
        if `which rate4site`.size==0
          raise "Rate4site not found"
        end
        alignfile = Tempfile.new('align')
        alignfile.puts self.output_clustal
        alignfile.flush
        outfile = Tempfile.new('rate4site')
        `rate4site -s #{alignfile.path} -o #{outfile.path} 2> /dev/null`
        output = outfile.grep(/^[^#]/).map{|line| line.chop}.delete_if{|i| i==''}
        output.map{ |line| 
          line.split[1]=~/B|X|Z/ ? nil : line.split[2].to_f
          }.compact
      end

      # seala
      # An interface to the SEALA software
      def seala(scorenum=1,suffix='bla')
        param2 = Tempfile.new 'param2'
        param2.puts "50 50\n1 1 1\n"
        param2.flush
        matrix = File.join(File.dirname(`which calcd`),"blosum62.#{suffix}")
        cons =
          "calcd -i #{alignfile.path} " + 
          "-w #{weightfile.path} " +
          "-f #{freqfile.path} -p #{param2.path} -m #{scorenum} " +
          "-t 2 -c #{matrix} -s #{self.keys.first}"
        output = `#{cons}`
        #puts output,'Press return' ; gets
        open('/tmp/cons','w'){|f| f.puts output}
        if output.nsplit.size<2
          raise "#{output}"
        end
        output.split("\n")[1..-1].map { |line|
          line=~/(B|X|Z)/ ? nil : line.split(",").last.to_f
        }.compact
      end

      # mihalek04seala
      # EThybrid - Evolutionary trace with shannon entropy
      def mihalek04seala
        exec = 
          "calci -i #{alignfile.path} " +
          "-p param3 -fp param4 -w 1 -value 2 -m 5 " + 
          "-t 2 -s #{self.keys.first}"
          #"-t 2 -s #{self.keys.first} -c blosum/blosum62.bla -e"
        ##################
        #puts exec
        #gets
        ##################
        output = `#{exec}`
        n = self.size
        output.split("\n").map { |line|
          if line=~/(B|X|Z)/ 
            nil 
          else 
            line.split("\t").last.to_f
          end
        }.compact.map{ |i|
          if i==100 then 'NaN'
          else i end } # (i-1).to_f / (n-1) end } # We normalize the score
      end

    end
  end
end

