function [vc,hw,wavefields] = mode_tracing(vc,hw,wavefields)
% mode tracing
   for freno=size(hw,2):-1:2
       I=1:size(hw,1);
       for id=1:find(~isnan(vc(:,freno-1)),1,'last')
           if isnan(real(hw(id,freno)))
              continue; 
           end
           cmin=min(abs(hw(:,freno-1)-hw(id,freno)));
           I(id)=find(abs(hw(:,freno-1)-hw(id,freno))==cmin,1);
       end

        index=find(diff(I)==0);
        I(index)=I(index)-1;

       vc(:,1:freno-1)=vc(I,1:freno-1);
       hw(:,1:freno-1)=hw(I,1:freno-1);
       
       for freno2=1:freno-1
           tempr=wavefields(freno2).ur;
           tempz=wavefields(freno2).uz;
           if size(tempr,2)<size(vc,2)
              tempr(:,size(vc,2))=0;
              tempz(:,size(vc,2))=0; 
           end
           wavefields(freno2).ur=tempr(:,I);
           wavefields(freno2).uz=tempz(:,I);
       end
       
   end
end

