require 'mathn'
require 'bigdecimal'
#require 'progressbar'

class Integer
  # calculates binomial coefficient of self choose k not recommended for large
  # numbers as binomial coefficients get large quickly... e.g. 100 choose 50 is
  # 100891344545564193334812497256
  def choose(k)
    return 0 if (k > self)
    n = self
    r = 1
    1.upto(k) do |d|
      r *= n
      r /= d
      n -= 1
    end
    return r
  end

  # Given n coin tosses, each with the probability p
  # for "success", what is the probability that
  # self is the number of successes?
  def binprob(n,p,normal=nil)
    #if n*p > 5 and n*(1-p) > 5
    if n*p > 10 and n*(1-p) > 10
    #if normal
      return (self.to_f).normalprob(n*p,n*p*(1-p))
    end
    unless p.is_a? Rational
      raise "Oh no you don't. Second argument must be a Rational." 
    end
    pnum = p.numerator.primes * self
    pden = p.denominator.primes * self
    p2num = (1-p).numerator.primes * (n-self)
    p2den = (1-p).denominator.primes * (n-self)
    numer = 
      (1..n).map{|i| i.primes}.flatten +
      pnum + p2num
    denom = 
      (1..self).map{|i| i.primes}.flatten +
      (1..n-self).map{|i| i.primes}.flatten +
      pden + p2den
    numer.delete(1)
    numer.sort!
    denom.delete(1)
    denom.sort!
    #puts numer.size,denom.size
    #p numer,denom
    arr = if numer.size <= denom.size
            numer.dup
          else
            denom.dup
          end
    arr.each do |i|
      if ni=numer.index(i) and di=denom.index(i)
        numer.delete_at(ni)
        denom.delete_at(di)
      end
    end
    #puts '*********'
    #puts numer.size,denom.size
    #p numer,denom
    #puts "old:\n"
    #puts n.choose(self) * p**self * (1-p)**(n-self)
    ans = (BigDecimal(numer.inject(1){|s,i| s=s*i}.to_s) /
           BigDecimal(denom.inject(1){|s,i| s=s*i}.to_s)).to_f
    if ans==0
      1e-300
    else
      ans
    end
  end

  def primes
    return [1] if self==1
    p = Prime.new
    n = self
    div = p.succ
    arr = []
    while n>1 do
      if (tmp = n / div).is_a? Fixnum
        n = tmp
        arr << div
      else
        div = p.succ
      end
    end
    arr
  end

end

class Float
  def ln
    Math::log(self)
  end

  def normalprob(m,d)
    ans = Math::exp(-(self-m)**2 / (2*d)) /
            (Math::sqrt(2*d*Math::PI))
    if ans == 0
      1e-300
    else
      ans
    end
  end
end

module Bio
  module Alignment
    class OriginalAlignment

      def seqs
        @seqs
      end

      def lockless99
        backgroundProbs = #by BLOSUM62
          {'A' => 74/1000,'R' => 52/1000,'N' => 45/1000,'D' => 54/1000,
           'C' => 25/1000,'Q' => 34/1000,'E' => 54/1000,'G' => 74/1000,
           'H' => 26/1000,'I' => 68/1000,'L' => 99/1000,'K' => 58/1000,
           'M' => 25/1000,'F' => 47/1000,'P' => 39/1000,'S' => 57/1000,
           'T' => 51/1000,'W' => 13/1000,'Y' => 32/1000,'V' => 73/1000}
        #backgroundProbs = #by Dayhoff (1978)
        #  {'A' => 86/1000,  'R' => 49/1000,  'N' => 43/1000,  'D' => 55/1000, 
        #   'C' => 29/1000,  'Q' => 39/1000,  'E' => 60/1000,  'G' => 84/1000, 
        #   'H' => 20/1000,  'I' => 45/1000,  'L' => 74/1000,  'K' => 66/1000, 
        #   'M' => 17/1000,  'F' => 36/1000,  'P' => 52/1000,  'S' => 70/1000, 
        #   'T' => 61/1000,  'W' => 13/1000,  'Y' => 34/1000,  'V' => 66/1000} 
        values = self.values.join
        n = self.size
        ncols = values.size / n
        raise "Not integer?" unless ncols.is_a? Fixnum
        #pbar = ProgressBar.new("Lockless99",self[self.keys.first].size)
        msaprob = backgroundProbs.inject({}) do |mp,(aa,bp)|
                    mp[aa] = (values.count(aa)/ncols).round.binprob(n,bp)
                    mp
                  end
        #puts "Average #{(values.count('A')/ncols).round}/#{n} times " +
        #     "(probability #{msaprob['A']})."
        ans =
          self.valind.map do |i|
            column = self.slice(i..i).values.join
            if (column.count('.').to_f / column.size) <= 0.5 #Maximum 50% gaps
              #pbar.set(i)
              Math::sqrt(backgroundProbs.inject(0) { |sqsum,(aa,bp)|
                           numer = column.count(aa).binprob(n, bp)
                           div = numer / msaprob[aa]
                           #if aa=='A'
                           #  puts "A found #{column.count('A')}/#{n} times " +
                           #       "(probability #{numer} => " +
                           #       "#{(Math::log(div))**2})"
                           #end
                           if numer==0
                             puts "Extraordinary event at #{i}:"
                             puts column
                             puts column.count(aa)
                             puts n
                             puts bp
                             puts msaprob[aa]
                             puts '*********'
                             gets
                           end
                           sqsum += (Math::log(div))**2 })
            else
              'nan'
            end
          end
        #pbar.finish
        ans
      end

    end
  end
end




