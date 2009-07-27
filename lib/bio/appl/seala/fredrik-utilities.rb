require 'rubygems'
require 'yaml'

def log(x)
  Math::log(x) / Math::log(10)
end

def ln(x)
  Math::log(x)
end

def log2(x)
  Math::log(x) / Math::log(2)
end

class Hash
  def hashmap
    self.inject({}) do |newhash, (k,v)|
      newhash[k] = yield(k, v)
      newhash
    end
  end
end

class String

  def fparse
    if self.split("\n").size > 1
      self
    else
      open(self).read
    end
  end

  def indices(substr)
    ans = []
    s = self
    startind = 0
    while i=s.index(substr)
      ans << i+startind
      s = s[i+1..-1]
      startind += i+1
    end
    ans
  end

  # Calculate percent identity to another string
  # by simply comparing each index's character
  # all the way until the end of self.
  def identity(other)
    n = self.split('').zip(other.split('')).sum{|i,j| i==j ? 1 : 0}
    n.to_f * 100 / [self.size, other.size].max
  end

end

class String
  def nsplit
    self.split(/[\n\t\s]*\n[\n\t\s]*/)
  end

  def to_r
    num,denom = self.split('/').map{|i| i.to_i}
    denom ||= 1
    Rational(num,denom)
  end
end

class Array

  def njoin
    self.join("\n")
  end
  def tjoin
    self.join("\t")
  end
  def sjoin
    self.join(' ')
  end

  def to_hash
    h = {}
    self.each_with_index do |si,i|
      h[i] = si
    end
    h
  end

  # Make it into a distribution with sum 1
  def distr!
    s = self.sum
    self.map! do |i|
      i / s
    end
  end

  def indices
    arr = []
    (0...self.size).each do |i|
      if yield self[i]
        arr << i
      end
    end
    arr
  end

end

module Enumerable

  def sum
    arr = if block_given?
            self.map{ |i| yield i }
          else
            self
          end
    if arr.first.is_a? Array
      arr.inject{ |s,i| s.zip(i).map{|j1,j2| j1+j2} }
    else
      arr.inject{ |s,i| s+i }
    end
  end

  def mean
    self.sum.to_f / self.size
  end

  def count(arg=nil)
    if block_given?
      self.sum{|i| yield(i) ? 1 : 0}
    else
      self.sum{|i| i==arg ? 1 : 0}
    end
  end

end


class String
  AAtranslation = {
    'A' => 'ALA',
    'R' => 'ARG',
    'N' => 'ASN',
    'D' => 'ASP',
    'C' => 'CYS',
    'E' => 'GLU',
    'Q' => 'GLN',
    'G' => 'GLY',
    'H' => 'HIS',
    'I' => 'ILE',
    'L' => 'LEU',
    'K' => 'LYS',
    'M' => 'MET',
    'F' => 'PHE',
    'P' => 'PRO',
    'S' => 'SER',
    'T' => 'THR',
    'W' => 'TRP',
    'Y' => 'TYR',
    'V' => 'VAL'}
  AAtranslation2 = {
    'ALA' => 'A',
    'ARG' => 'R',
    'ASN' => 'N',
    'ASP' => 'D',
    'CYS' => 'C',
    'GLU' => 'E',
    'GLN' => 'Q',
    'GLY' => 'G',
    'HIS' => 'H',
    'ILE' => 'I',
    'LEU' => 'L',
    'LYS' => 'K',
    'MET' => 'M',
    'PHE' => 'F',
    'PRO' => 'P',
    'SER' => 'S',
    'THR' => 'T',
    'TRP' => 'W',
    'TYR' => 'Y',
    'VAL' => 'V'}
  def threeLetterAA()
    if ans = AAtranslation[self.upcase]
      ans
    else
      raise ArgumentError, "Not a one letter amino acid"
    end
  end

  def oneLetterAA()
    if ans = AAtranslation2[self.upcase]
      ans
    else
      raise ArgumentError, "Not a three letter amino acid"
    end
  end
end
