require 'rsruby'
require 'tempfile'

STDERR.puts "   (Starting R interpreter...)"
R = RSRuby.instance # keep a constant R interpreter
R.matrix.autoconvert(RSRuby::NO_CONVERSION)
#R.as_dist.autoconvert(RSRuby::NO_CONVERSION)

#Dirty hack of getting some R code:
f = Tempfile.new('neumannR')
f.puts 'neumann <- function(obs,densityfile,logbase=0) {
    densitymatrix <- as.matrix(read.table(densityfile))
    w <- diag(obs) %*% densitymatrix
    w <- w / sum(diag(w))
    e <- eigen(w, only.values=TRUE)$values
    if (logbase>0)
      { factor = log(logbase)
        - sum(e * log(e) / factor,na.rm=TRUE) }
    else
        - sum(e * log(e),na.rm=TRUE)
}
'
f.flush
R.source f.path
R.library 'cluster'

module L

  # --- stats
  def L.mean items
    R.mean items
  end
  
  def L.stddev items
    R.sd items
  end
  
  def L.saveplot(filename, data, params=nil)
    if File.exists? filename
      raise ArgumentError, "File already exists"
    end
    R.pdf filename 
    R.plot(data, params)
    R.eval_R("dev.off()")
  end

  def L.rank data
    R.rank(data).map do |i|
      i.to_i
    end
  end

  def L.neumann(profile,densityfile,logbase=0)
    arr = Bio::Alignment::Alphabet.map { |a|
            profile[a] ? profile[a] : 0 }
    R.neumann(arr,densityfile,logbase)
  end

end

class Matrix

  def size
    [self.row_size, self.column_size]
  end

  def flatten
    self.to_a.flatten
  end

  def to_R
    R.matrix(self.flatten.map{|i|i.to_f},
             :nrow=>self.row_size,:ncol=>self.row_size)
  end

  def eigenfunc
    raise ArgumentError, "Need a square matrix" unless self.square?
    m = self.to_R
    R.eigen(m)
  end

  def eigen
    self.eigenfunc['vectors'].map do |i|
      Matrix.columns([i])
    end
  end

  def eigenval
    self.eigenfunc['values']
  end

  def to_s
    (0...self.row_size).map do |i|
      self.row(i).to_a.map{|j| j.to_s.rjust(5)}.join(" ")
    end.join("\n")
  end
end

class Array
  def rank
    L.rank self
  end
end
