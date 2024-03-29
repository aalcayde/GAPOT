classdef PowerTimeDomain
   properties
      Va
      Vb
      Vc
      gaV2a
      gaV2b
      gaV2c
      gaV2
      Varms
      Vbrms
      Vcrms
      Vrms
      gaVa
      gaVb
      gaVc
      gaV
      gaI
      gaHV
      gaVrms
      HVa
      HVb
      HVc
      gaMa
      gaMb
      gaMc
      gaM
      Mpa
      Mpa_osci
      Mpb
      Mpc
      gaMp
      Mqa
      Mqa_osci
      Mqb
      Mqc
      gaMq
      mMpa
      mMqa
      P
      Q
      Ia
      Ib
      Ic
      gaIp
      gaIq
      gaIf
      gaIx
      gaIbr
      gaIu
      
      gaIa
      gaIb
      gaIc
      Ipa
      Ipb
      Ipc
      Iqa
      Iqb
      Iqc
      Ifa
      Ifb
      Ifc
      Ixa
      Ixb
      Ixc
      Ibra
      Ibrb
      Ibrc
      Iua
      Iub
      Iuc
      Iarms
      Ibrms
      Icrms
      Irms
      Iparms
      Iqarms
      Ifarms
      Ixarms
      Ibrarms
      HIa
      HIb
      HIc
      SinglePhase
      LoadedDataSinglePhase
      LoadedData
      SimulationTime
      W
      Fs
   end
   methods
      function r = roundOff(obj)
         r = round([obj.Value],2);
      end
      function r = multiplyBy(obj,n)
         r = [obj.Value] * n;
      end
   end
end