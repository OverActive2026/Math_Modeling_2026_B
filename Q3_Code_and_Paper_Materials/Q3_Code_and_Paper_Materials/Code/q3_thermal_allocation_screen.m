function screening=q3_thermal_allocation_screen()
% Generate preheat power proposals from the baseline mean thermal network.
% This is an optimizer screening model only: no screening row is a Q3 result.
% Water/ice latent heat remains in all final physical candidate evaluations.
    area=25e-4;
    delta=[150e-6;3.4e-6;12e-6;11.3e-6;150e-6];
    rho=[185;970;2150;970;185]; cp=[545;240;1050;240;545];
    k=[.3;.27;.24;.27;.3];
    meaCapacity=sum(rho.*cp.*delta);
    scale=[.75;.5;.5;.5;.75];
    C=[area*(meaCapacity+1980*766*.004*scale);98.75;98.75];
    resistance=sum(delta./k);
    Gcc=area/(resistance+.002/95);
    Gend=area/(.5*resistance+.002/95+.005/15);
    Genv=area/(.005/15+1/40);
    K=zeros(7);
    for j=1:4, add_edge(j,j+1,Gcc); end
    add_edge(1,6,Gend); add_edge(5,7,Gend);
    K(6,6)=K(6,6)+Genv; K(7,7)=K(7,7)+Genv;
    A=-diag(1./C)*K;
    B=diag(1./C)*[25*eye(5);zeros(2,5)];
    block=[A B;zeros(5,12)];
    times=[25:.05:50,50.5:.5:180].';
    E=nan(size(times)); powers=nan(numel(times),5);
    options=optimoptions('linprog','Display','none');
    for j=1:numel(times)
        response=expm(block*times(j)); H=response(1:5,8:12);
        [q,fval,flag]=linprog(ones(5,1),-H,-30.01*ones(5,1), ...
            [],[],zeros(5,1),ones(5,1),options);
        if flag>0
            powers(j,:)=q.'; E(j)=25*times(j)*fval;
        end
    end
    screening=table(times,powers,E,'VariableNames',{'timeS','q','energyEstimateJ'});
    outDir=fullfile(fileparts(mfilename('fullpath')),'results');
    model_id='thermal_allocation_screen_only';
    [~,idx]=min(E); best=screening(idx,:);
    save(fullfile(outDir,'thermal_allocation_screen.mat'), ...
        'screening','best','model_id','C','Gcc','Gend','Genv');
    disp(best);
    function add_edge(i,j,g)
        K(i,i)=K(i,i)+g; K(j,j)=K(j,j)+g;
        K(i,j)=K(i,j)-g; K(j,i)=K(j,i)-g;
    end
end
