
module Bio
  module Alignment
    class OriginalAlignment
      
      Taylorsets =
        %w(RKH
        DERKH
        DERK
        DE
        SNCEQ
        SNDEQR
        SNDEQRKH
        SNDEQRKHP
        PTSNDEQRKH
        TSNDEQRKHWY
      
        PTSNDEQRKHWY
        TSNDEQRKH
        PTSNDEQ
        TSNDEQ
        PTSND
        TSND
        SND
        AGS
        AGTSND
        PAGTSND
      
        AGTSNDEQ
        PAGTSNDEQ
        AGTSNDEQRK
        PAGTSNDEQRK
        AGTSNDEQRKHWY
        VCAGTSNDEQRKHWY
        PVCAGTSNDEQRKHWY
        VCAGTSNDEQRK
        PVCAGTSNDEQRK
        VCAGTSNDEQR
      
        PVCAGTSND
        VCAGTSND
        VCAGTS
        VCAGT
        VCAG
        VCAGP
        VCAGPT
        LIVCAGPT
        LIVCAGP
        LIVCAG
      
        LIV
        FMLIV
        FMLIVCAG
        PFMLIVCAG
        FMLIVCAGT
        PFMLIVCAGT
        FMLIVCAGTK
        FMLIVCAGTKSND
        FMLIVCAGTKSNDP
        HWYFMLIVCAGTKSNDP
      
        HWYFMLIVCAGTKSND
        HWYFMLIVCAGTK
        HWYFMLIVCAGTKP
        HWYFMLIVCAG
        HWYFMLIV
        HWYFM
        HWYF
        RKHWYF
        QERKHWY
        QERKHWYFM
        QERKHWYFMLI
        )

      def taylor86
        self.valind.map do |ci|
          obs = self.slice(ci..ci).values.uniq.
                reject{ |i| !Alphabet.member? i}.join
          memberof = Taylorsets.find_all{|i| i =~ /([#{obs}].*){#{obs.size}}/}
          if memberof.size>0
            memberof.map{|i| i.size}.min
          else
            20
          end
        end
      end

    end
  end
end
