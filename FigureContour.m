function h=FigureContour(z,P,order,tresh_z,Type,SaveFig,MainTitle,Tag);
    scrsz = get(0,'ScreenSize');
    ox=order(1); oy=order(2); or=order(3); oc=order(4);
    nx=P(ox).n; ny=P(oy).n; nr=P(or).n; nc=P(oc).n;
    x=P(ox).a; y=P(oy).a; r=P(or).a; c=P(oc).a;
    lx=P(ox).l; ly=P(oy).l; lr=P(or).l; lc=P(oc).l;
    ux=P(ox).u; uy=P(oy).u; ur=P(or).u; uc=P(oc).u;

    figure('Name',Type,'Position',[0 0 scrsz(3) scrsz(4)]);
    levels=[1E-7 1E-6 1E-5 1E-4 1E-3 1E-2 3E-2 6E-2 1E-1 3E-1 6E-1 1E0 1E1 1E2 1E3 1E4 1E5 1E6 1E7 1E8 1E9 1E10 1E11 1E12];
    subplot1(nr,nc, 'XTickL', 'Margin', 'YTickL', 'Margin','YScale','linear');
    for ir=1:nr,
        for ic=1:nc,
            subplot1(ic+nc*(ir-1));
            zvalues=squeeze(z(:,:,ir,ic)');
            contour(x,y,zvalues,levels,'b','ShowText','on'); hold on;
            %set(gca,'YScale','log')
            contour(x,y,zvalues,[tresh_z 1E30],'k','LineWidth',2);
            %contourf(x,y,zvalues,[-1E30 1E-3],'r','LineWidth',2);
            if(min(x)<0.5*max(x)), minx=0; else minx=min(x); end
            if(min(y)<0.5*max(y)), miny=0; else miny=min(y); end
            xlim([minx max(x)]); ylim([miny max(y)]);
            if(ir==nr), xlabel([lx ' (' ux ')']); end
            if(ic==1), ylabel({[lr ' = ' P(or).s{ir} ' ' ur];[ly ' (' uy ')']}); end
            if(ir==1), title([lc ' = ' P(oc).s{ic} ' ' uc]); end
            if(strcmp(Type,'AMBIENT')), text(max(x)/2,max(y)/2,num2str(max(max(zvalues)),3)); end
            %pcolor(x,y,zvalues); shading interp; colormap(pink);
            hold off;
        end
    end
    h=0;
    suptitle([Type ' - ' MainTitle]);
    NameFig=['Results\' Type '_' P(ox).l '_' P(oy).l '_' P(or).l '_' P(oc).l '_' Tag];
    if SaveFig, save_figure(NameFig); end
end