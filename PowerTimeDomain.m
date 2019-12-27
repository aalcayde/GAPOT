classdef PowerTimeDomain
   properties
      Va
      Vb
      Vc
      gaV2a
      gaV2b
      gaV2c
      Varms
      Vbrms
      Vcrms
      gaVa
      gaVb
      gaVc
      HVa
      HVb
      HVc
      gaMa
      gaMb
      gaMc
      Mpa
      Mpb
      Mpc
      Mqa
      Mqb
      Mqc
      mMpa
      mMqa
      P
      Q
      Ia
      Ib
      Ic
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