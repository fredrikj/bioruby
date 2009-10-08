module Bio
  module Alignment
    class OriginalAlignment

      BLOSUM62freq =
        {'A' => 74/1000,
         'R' => 52/1000,
         'N' => 45/1000,
         'D' => 54/1000,
         'C' => 25/1000,
         'Q' => 34/1000,
         'E' => 54/1000,
         'G' => 74/1000,
         'H' => 26/1000,
         'I' => 68/1000,
         'L' => 99/1000,
         'K' => 58/1000,
         'M' => 25/1000,
         'F' => 47/1000,
         'P' => 39/1000,
         'S' => 57/1000,
         'T' => 51/1000,
         'W' => 13/1000,
         'Y' => 32/1000,
         'V' => 73/1000}

      def pdist
        return @pdist if @pdist
        str = self.values.join.split('').reject{|i| !Alphabet.member? i}
        strsize = str.size
        @pdist = {}
        Alphabet.each do |a|
          @pdist[a] = Rational(str.count(a), strsize)
        end
        @pdist
      end

      def mihalek04
        self.evoltrace 'shannon'
      end

      def zhang08
        self.evoltrace 'neumann'
      end

      def shannon
        self.valind.map do |ci|
          entropy(self.profile[ci], Alphabet,20)
        end
      end

      def shannonw
        self.setweights! 3
        ans = 
          self.valind.map do |ci|
            entropy(self.profile[ci], Alphabet,20)
          end
        self.setweights! 1
        ans
      end

      def neumann
        self.valind.map do |ci|
          if (p=self.profile[ci]).size<1
            0.0
          else
            matrix = File.join(File.dirname(`which calcd`),"blosum62.qij2")
            arr = Alphabet.map { |a| p[a] ? p[a] : 0 }
            R.neumann(arr,matrix,20)
          end
        end
      end

      def wang06(backdist=BLOSUM62freq)
        self.valind.map do |ci|
          relentropy(self.profile[ci], backdist,20)
        end
      end

      def wang06w(backdist=BLOSUM62freq)
        self.setweights! 3
        ans =
          self.valind.map do |ci|
            relentropy(self.profile[ci], backdist,20)
          end
        self.setweights! 1
        ans
      end

      def capra07(backdist=BLOSUM62freq)
        self.valind.map do |ci|
          jsd(self.profile[ci], backdist,20)
        end
      end

      def capra07w(backdist=BLOSUM62freq)
        self.setweights! 3
        ans = 
          self.valind.map do |ci|
            jsd(self.profile[ci], backdist,20)
          end
        self.setweights! 1
        ans
      end

      def mihalek07(backdist=nil)
        backdist ||= readMihalek
        self.valind.map do |ci|
          jointrelentropy(self.slice(ci..ci).values, backdist)
        end
      end

    end
  end
end

def to_freq(observations,alphabet)
  if observations.is_a? Hash
    observations
  elsif observations.class.to_s =~ /(String)|(Array)/
    observations = observations.split('') if observations.is_a? String
    alphabet = alphabet.split '' if alphabet.is_a? String
    n = observations.reject{ |i| !alphabet.member? i}.size
    alphabet.inject Hash.new do |h,i| 
      h[i] = Rational(observations.count(i), n) ; h
    end
  end
end

def relentropy(observations, backdist, logbase=nil)
  alphabet = backdist.keys
  freq = to_freq(observations, alphabet)
  if logbase
    factor = ln(logbase)
    alphabet.sum{ |i|
      (!freq[i] or freq[i]==0) ? 0 : freq[i]*ln(freq[i]/backdist[i])/factor }
  else
    alphabet.sum{ |i|
      (!freq[i] or freq[i]==0) ? 0 : freq[i] * ln(freq[i] / backdist[i]) }
  end
end

def entropy(observations,alphabet,logbase=nil)
  alphabet = alphabet.split '' if alphabet.is_a? String
  alphabet = alphabet.keys if alphabet.is_a? Hash
  freq = to_freq(observations, alphabet)
  if logbase
    factor = ln(logbase)
    alphabet.sum{ |i| 
      (!freq[i] or freq[i]==0) ? 0 : - freq[i] * ln(freq[i]) / factor }
  else
    alphabet.sum{ |i| 
      (!freq[i] or freq[i]==0) ? 0 : - freq[i] * ln(freq[i]) }
  end
end

def jsd(observations, backdist,logbase=nil)
  p = to_freq(observations,backdist.keys)
  avg = backdist.hashmap { |aa, val| [val, p[aa]].mean }
  [ relentropy(p,avg,logbase), relentropy(backdist,avg,logbase) ].mean
end

def jointrelentropy(obs,backdist)
  obs = obs.split('') if obs.is_a? String
  alphabet = Bio::Alignment::Alphabet
  n = obs.reject{|i| !alphabet.member? i}.size
  return 'nan' if n<2
  maxpairs = n * (n-1) / 2
  dist = {}
  alphabet.each_with_index do |a1,i|
    alphabet[0..i].each do |a2|
      num = (a1==a2 ? obs.count(a1).choose(2) : obs.count(a1)*obs.count(a2))
      dist[a1+a2] = Rational(num, maxpairs)
    end
  end
  relentropy(dist,backdist)
end

def readMihalek
  str = 
'A R N D C Q E G H I L K M F P S T W Y V
0.1681
0.0451 0.1917
0.0402 0.0376 0.1401
0.0467 0.0455 0.0704 0.2289
0.0222 0.0208 0.0282 0.0178 0.3769
0.0340 0.0530 0.0469 0.0326 0.0195 0.1177
0.0527 0.0521 0.0516 0.0813 0.0529 0.0940 0.1696
0.0848 0.0536 0.0931 0.0576 0.0444 0.0600 0.0476 0.2842
0.0201 0.0253 0.0392 0.0254 0.0333 0.0257 0.0232 0.0201 0.2344
0.0466 0.0357 0.0493 0.0324 0.0293 0.0369 0.0303 0.0263 0.0352 0.1416
0.0626 0.0632 0.0565 0.0593 0.0594 0.0622 0.0069 0.0011 0.0640 0.1208 0.2021
0.0585 0.0857 0.0636 0.0564 0.0458 0.0650 0.0738 0.0351 0.0441 0.0441 0.0195 0.1443
0.0221 0.0241 0.0275 0.0130 0.0180 0.0310 0.0151 0.0178 0.0731 0.0273 0.0345 0.0234 0.0612
0.0317 0.0375 0.0368 0.0253 0.0261 0.0325 0.0260 0.0250 0.0454 0.0500 0.0586 0.0361 0.0596 0.1903
0.0359 0.0374 0.0350 0.0398 0.0250 0.0467 0.0461 0.0372 0.0552 0.0388 0.0288 0.0464 0.0289 0.0251 0.2932
0.0719 0.0547 0.0587 0.0527 0.0377 0.0671 0.0542 0.0385 0.0554 0.0462 0.0037 0.0573 0.0490 0.0481 0.0450 0.1323
0.0514 0.0375 0.0489 0.0401 0.0453 0.0458 0.0450 0.0348 0.0458 0.0452 0.0457 0.0495 0.0503 0.0417 0.0509 0.0677 0.1167
0.0106 0.0135 0.0079 0.0079 0.0299 0.0394 0.0109 0.0101 0.0426 0.0089 0.0121 0.0123 0.3051 0.0463 0.0184 0.0097 0.0348 0.1789
0.0313 0.0336 0.0274 0.0240 0.0267 0.0514 0.0274 0.0273 0.0592 0.0378 0.0379 0.0300 0.0477 0.0901 0.0255 0.0312 0.0281 0.1588 0.1440
0.0633 0.0524 0.0411 0.0430 0.0409 0.0388 0.0396 0.0015 0.0333 0.1173 0.0012 0.0095 0.0714 0.0678 0.0407 0.0189 0.0746 0.0420 0.0608 0.1419'
  #str = open('Mihalek07-Q').read
  keys = str.split("\n").first.split
  (0...keys.size).zip(keys).inject({}) do |hash,(row,k)|
    (0..row).zip(keys).each do |col,k2|
      hash[k+k2] = str.split("\n")[row+1].split[col].to_f
    end
    hash
  end
end
