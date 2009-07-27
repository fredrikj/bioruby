require 'set'

class Rtree < Array
  #def initialize(arg)
  #  super arg.split("\n").map{|i| i.split.map{|k| k.to_i}}
  #end

  def children(node)
    if (node < 0)
      raise ArgumentError, "No such element" if -node>self.size+1
      Set[-node]
    else
      raise ArgumentError, "No such node" if node>self.size
      fg = self[node-1]
      self.children(fg[0]) + self.children(fg[1])
    end
  end
      
  # Remember - the nodes are numbered from bottom to up
  def groups(node)
    if node>self.size or node<1
      raise ArgumentError, "No such node"
    end
    g = Set.new((1..self.size+1))
    k = node
    while k>0 do
      c = self.children(k)
      if !(g.find{|ki| ki.is_a?(Set) and c.subset?(ki)})
        g -= c
        g += [c]
      end
      k -= 1
    end
    g
  end
end
