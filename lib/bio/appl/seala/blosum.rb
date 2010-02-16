require 'matrix'

# a = Blosum.new string
# a['A']['B'] => 1

class Blosum < Hash
  attr_reader :keys
  def initialize(arg)
    @keys, *values = arg.lines.grep(/^[^#]/).map{|i| i.split}
    @ndec = (values[0][1]=~/\.(\d*)/) ? $1.size : 0
    @keys.each_with_index do |key,i|
      self[key] = 
        @keys.zip((0...@keys.size).to_a).inject({}) { |hash,(key2,j)|
          row = [i,j].max
          col = [i,j].min
          hash[key2] = values[row][col].to_f; hash }
    end
  end

  def to_s
    (['   ' + @keys.join('  ') ] +
      @keys.map{ |key| 
        key + ' ' + @keys.map{ |key2| 
                      ("%.#{@ndec}f"%self[key][key2]).rjust(2) 
                    }.join(' ') 
      }).join("\n")
  end

  def karlin!
    diag = @keys.inject({}){|hash,key| hash[key] = self[key][key].abs; hash}
    @keys.each do |k1|
      @keys.each do |k2|
        self[k1][k2] = self[k1][k2] / Math::sqrt(diag[k1] * diag[k2])
      end
    end
    self
  end

  def translate!(arg)
    @keys.each do |k1|
      @keys.each do |k2|
        self[k1][k2] += arg
      end
    end
    self
  end
end


class Blosum2 < Hash
  attr_reader :keys
  attr_accessor :ndec
  def initialize(arg)
    # First, remove all comments
    lines = arg.split("\n").reject{|line| line=~/#/ or line.size==0}
    # First line is expected to list all keys
    @keys = lines.shift.split
    raise "Ill formed matrix" if @keys.size != lines.size
    # Now, parse all lines
    values = lines.map { |i| i.split }
    (0...lines.size).each do |i|
      if values[i][0] =~ /[A-Z]/
        key = values[i].shift
        raise "Bad key value (#{key}) on line #{i+1}" if key!=@keys[i]
      end
    end
    @ndec = (values[0][0]=~/\.(\d*)/) ? $1.size : 0
    @keys.each_with_index do |key,i|
      @keys.each_with_index do |key2,j|
        row = [i,j].max
        col = [i,j].min
        self[key,key2] = values[row][col].to_f
        self[key2,key] = values[row][col].to_f
      end
    end
  end

  def []=(arg1,arg2,val)
    super(arg1+arg2,val)
  end

  def [](arg1,arg2=nil)
    if arg2
      super(arg1+arg2)
    else
      @keys.inject({}) do |hash,k|
        hash[k] = self[arg1,k]
        hash
      end
    end
  end

  def to_s
    (['   ' + @keys.join('  '+' '*@ndec) ] +
      @keys.map{ |key| 
        key + ' ' + @keys.map{ |key2| 
                      ("%.#{@ndec}f"%self[key,key2]).rjust(2) 
                    }.join(' ') 
      }).join("\n")
  end

  def karlin!
    diag = @keys.inject({}){|hash,key| hash[key] = self[key,key].abs; hash}
    @keys.each do |k1|
      @keys.each do |k2|
        self[k1,k2] = self[k1,k2] / Math::sqrt(diag[k1] * diag[k2])
      end
    end
    self
  end

  def translate!(arg)
    @keys.each do |k1|
      @keys.each do |k2|
        self[k1,k2] += arg
      end
    end
    self
  end
end


