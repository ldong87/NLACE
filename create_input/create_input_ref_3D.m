% Loading the grid coordinates and the measured displacement at the right format (cf README)
%filename='CC024summedDisplacements21_35f.mat';

%clc
close all
clear all

filename='test3D-new.mat';% read xDisp, yDisp, zDisp, pDisp, xCoor, yCoor, zCoor, nDim, comment
filename='dseries_accum_110207_1f.mat';% read xDisp, yDisp, zDisp, pDisp, xCoor, yCoor, zCoor, nDim, comment
filename='nonlinph-accum-20110218-ud3f.mat';% read xDisp, yDisp, zDisp, pDisp, xCoor, yCoor, zCoor, nDim, comment
filename='nonlinph-accum-20110302-ud3f-2.mat';% read xDisp, yDisp, zDisp, pDisp, xCoor, yCoor, zCoor, nDim, comment
filename='nonlinph-accumfilt-20110510.ud3f.mat';% read xDisp, yDisp, zDisp, pDisp, xCoor, yCoor, zCoor, nDim, comment
filename='nonlinph-accumfilt-20110731.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a20d10-20110808.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a100d10-20110811.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a1000d10-20110811.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a100000d20-20110824.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a1ep7d20-20110825.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a1ep7d1000-20110825.ud3f.mat';% incremental data is filtered and them accumulated
filename='nonlinph-accumfilt_a1ep7d20-20110825.ud3fpfall.mat';% filtered data by svd
filename='nonlinph-accumfilt_a1ep7d20-20110825.ud3fpfSmallStrain.mat';% fileter data by svd, small strain only

filename='/bigtmp/jfdord/Backup/london_data/inversions/incrementalDisp-slip1geoCylDispCart.mat';% structured irregular grid with cartesian displacements given at each node
%filename='/bigtmp/jfdord/Backup/london_data/inversions/incrementalDisp-slip5_1geoCylDispCart.mat';% structured irregular grid with cartesian displacements given at each node
filename='/bigtmp/jfdord/Backup/london_data/inversions/incrementalDisp-slip1.mat';% structured irregular grid with cartesian displacements given at each node
%filename='/bigtmp/jfdord/Backup/london_data/inversions/incrementalDisp-slip5_1.mat';% structured irregular grid with cartesian displacements given at each node
filename='/bigtmp/jfdord/Backup/london_data/inversions/incrementalDisp-slip1geoCylDispCartOrtho.mat';% structured irregular grid with cartesian displacements given at each node
%filename='/bigtmp/jfdord/Backup/london_data/inversions/incrementalDisp-slip5_1geoCylDispCartOrtho.mat';% structured irregular grid with cartesian displacements given at each node

disp(['filename is : ' filename]);

% select a series of loading steps
Ln=[1];%[10]% frame 1 for shear modulus
% select every "reduce" measured point
reduce=4;%3%2%18 for tiny size
% if the data on its boundary is inaccurate, define margins (in number fo lines/columns)
topMargin=26;%9;%2;% closest to the ultrasound transducer
bottomMargin=36;%20;%1;% far away from transducer
%sideMargin=1;% left and right sides are treated equally
elevMargin=20;%10;%1;
lateMargin=1;
% choose a type of element
elType=305;%32-> filter to smooth lateral+elev while respecting incompressibility;%305->3D incompressible;

addElem611=false;
notop = true ; %this will add/remove then top springs
notopandbottom = false;% this will add/remove the top and bottom springs

addElem307=false;
surfaceForce=50*50;%mm2
totalForce=647.0;%mN %frame20
%totalForce=1456.0;%mN %frame25

geoCylDispCart=true;% the displacements are still given in the cartesian xyz coordinates; the grid is structured but irregular
enforceOrthogonality=true;

%plotBCs=false;% for debug or information

%----------------------- creating the input file ---------------%
%if plotBCs  disp('check the filename hereafter...'); load target1_compression_accumf.mat; xDispi=xDisp; end
if ((elType==31)||(elType==32))&&(size(Ln,2)~=1)
  disp(['when using the filter(elType==31), only one displacement can be given']);
  return
end
load(filename);

%if plotBCs  xDispNlm=xDisp; end

% size of the data
ex = size(xDisp,2);
ey = size(xDisp,1);
ez = size(xDisp,3);
% apply the margins
%ex=ex-sideMargin;
ex=ex-lateMargin;
ey=ey-bottomMargin;
%ez=ez-sideMargin;
ez=ez-elevMargin;
%initial index for choosing data on x- and y-axes
%sx = 1+sideMargin;
sx = 1+lateMargin;
sy = 1+topMargin;
%sz = 1+sideMargin;
sz = 1+elevMargin;
% vector of the selected indexes
ax=sx:reduce:ex;
ay=sy:reduce:ey;
az=sz:reduce:ez;

% Normalizing measurement with its absolute minimum (useful for linear reconstruction)
%for i=1:1:size(Ln,2)
%  factor = max(max(abs(yDisp(:,:,Ln)),[],1),[],2);    
%  xDisp(:,:,Ln)=xDisp(:,:,Ln)/factor;
%  xDisp(:,:,Ln)=yDisp(:,:,Ln)/factor;
%end


if false
  % plot the displacements along some lines in the domain and, if chosen, smooth some measured displacements

  % select all indexes for ploting displacements
  bor=3;% border to take a few more nodes and perform a good quality smoothing (afterwards)
  ax2=max([ax(1)-bor,1],[],2):1:min([ax(end)+bor,size(xDisp,2)],[],2);
  ay2=max([ay(1)-bor,1],[],2):1:min([ay(end)+bor,size(xDisp,1)],[],2);
  bxmin=ax(1)-ax2(1);% keep the difference in mind to plot and save the right data
  bxmax=ax2(end)-ax(end);
  bymin=ay(1)-ay2(1);
  bymax=ay2(end)-ay(end);
  xDisp_mod=xDisp(ay2,ax2,Ln);
  yDisp_mod=yDisp(ay2,ax2,Ln);
  xCoor_mod=xCoor(ay2,ax2);
  yCoor_mod=yCoor(ay2,ax2);
  m=size(xDisp_mod,1);
  n=size(xDisp_mod,2);
  % check the noise on the x- and y-component of the displacement on the boundaries of the domain
  for lln=1:1:size(Ln,2)
    figure
    set(gcf,'color','w');
    subplot(1,2,1)
    hold on
    plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:),'k--');
    plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:)+xDisp_mod(1+bymin,:,lln),'b-');
    plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:),'k--');
    plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:)+xDisp_mod(floor(m/2),:,lln),'b-');
    plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:),'k--');
    plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:)+xDisp_mod(m-bymax,:,lln),'b-');
    plot(xCoor_mod(:,1+bxmin),yCoor_mod(:,1+bxmin),'k--');
    plot(xCoor_mod(:,1+bxmin)+xDisp_mod(:,1+bxmin,lln),yCoor_mod(:,1+bxmin),'r-');
    plot(xCoor_mod(:,floor(n/2)),yCoor_mod(:,floor(n/2)),'k--');
    plot(xCoor_mod(:,floor(n/2))+xDisp_mod(:,floor(n/2),lln),yCoor_mod(:,floor(n/2)),'r-');
    plot(xCoor_mod(:,n-bxmax),yCoor_mod(:,n-bxmax),'k--');
    plot(xCoor_mod(:,n-bxmax)+xDisp_mod(:,n-bxmax,lln),yCoor_mod(:,n-bxmax),'r-');
    title('Initial data: lateral disp. along some lines of the domain');
    legend('lines in the domain');
    axis equal
    subplot(1,2,2)
    hold on
    plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:),'k--');
    plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:)+yDisp_mod(1+bymin,:,lln),'b-');
    plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:),'k--');
    plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:)+yDisp_mod(floor(m/2),:,lln),'b-');
    plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:),'k--');
    plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:)+yDisp_mod(m-bymax,:,lln),'b-');
    plot(xCoor_mod(:,1+bxmin),yCoor_mod(:,1+bxmin),'k--');
    plot(xCoor_mod(:,1+bxmin)+yDisp_mod(:,1+bxmin,lln),yCoor_mod(:,1+bxmin),'r-');
    plot(xCoor_mod(:,floor(n/2)),yCoor_mod(:,floor(n/2)),'k--');
    plot(xCoor_mod(:,floor(n/2))+yDisp_mod(:,floor(n/2),lln),yCoor_mod(:,floor(n/2)),'r-');
    plot(xCoor_mod(:,n-bxmax),yCoor_mod(:,n-bxmax),'k--');
    plot(xCoor_mod(:,n-bxmax)+yDisp_mod(:,n-bxmax,lln),yCoor_mod(:,n-bxmax),'r-');
    title('Axial disp. along some lines of the domain (if no disp., the dashed and colored lines superpose)');
    legend('Lines in the domain');
    axis equal;
 %   figure
 %   subplot(1,2,1)
 %   surf(xCoor_mod,yCoor_mod,xDisp_mod(:,:,lln),'edgealpha',0);
 %   view(0,90);
 %   axis equal;
 %   title('Initial data: lateral displacement');
 %   subplot(1,2,2)
 %   surf(xCoor_mod,yCoor_mod,yDisp_mod(:,:,lln),'edgealpha',0);
 %   view(0,90);
 %   axis equal;
 %   title('axial displacement');
  end
  if false% use filters
    if false% laplace smoothing
      disp('Laplace smoothing everywhere');
      for i=1:1:30% smooth several times because it is laplace
        for lln=1:1:size(Ln,2)
          disp_tmp=zeros(m+2,n+2);
          disp_tmp(2:(m+1),2:(n+1))=xDisp_mod(:,:,lln);
          disp_tmp(2:(m+1),1)=xDisp_mod(:,1,lln);
          disp_tmp(2:(m+1),n+2)=xDisp_mod(:,n,lln);
          disp_tmp(1,:)=disp_tmp(2,:);
          disp_tmp(m+2,:)=disp_tmp(m+1,:);
          xDisp_mod(:,:,lln)=(disp_tmp(2:(m+1),2:(n+1))+...
                              disp_tmp(1:m,2:(n+1))+disp_tmp(3:(m+2),2:(n+1))+...
	    		      disp_tmp(2:(m+1),1:n)+disp_tmp(2:(m+1),3:(n+2)))/5;
        end
      end
    elseif false% laplace smoothing on the boundaries of the domain
      disp('Laplace smoothing on the boundaries');
      nLines=20;% number of line where the laplace filter is applied
      nTrans=4;% number of line for transitioning from 1 to 0
      filtre=ones(m,n);
      for i=1:1:nTrans
        filtre((nLines+i):(m-nLines-i),(nLines+i):(n-nLines-i))=(1-i/nTrans);
      end
      for i=1:1:30% smooth several times because it is laplace
        for lln=1:1:size(Ln,2)
          disp_tmp=zeros(m+2,n+2);
          disp_tmp(2:(m+1),2:(n+1))=xDisp_mod(:,:,lln);
          disp_tmp(2:(m+1),1)=xDisp_mod(:,1,lln);
          disp_tmp(2:(m+1),n+2)=xDisp_mod(:,n,lln);
          disp_tmp(1,:)=disp_tmp(2,:);
          disp_tmp(m+2,:)=disp_tmp(m+1,:);
          xDisp_mod(:,:,lln)=(disp_tmp(2:(m+1),2:(n+1))+filtre.*...
                              (disp_tmp(1:m,2:(n+1))+disp_tmp(3:(m+2),2:(n+1))+...
	    		       disp_tmp(2:(m+1),1:n)+disp_tmp(2:(m+1),3:(n+2))))./(1+4*filtre);
        end
      end
    elseif false% fourier filter
      % articially augment the size of the data, make it continuous to avoid jumps and troublesome ffts
      nn=7;% add nn ficticious cells all around the domain
      mm=2;% among the nn ficticious cells, mm will have the values at the boundary of the domain (mm<nn)
           %  the others create a linear transition
      nKeep=10;% number of most energetic modes to keep in every direction (all others are set to 0)
      for lln=1:1:size(Ln,2)
        disp_tmp=zeros(m+2*nn,n+2*nn);
        disp_tmp((1+nn):(m+nn),(1+nn):(n+nn))=xDisp_mod(:,:,lln);
        for j=1:1:nn% mm lines have the same value as the boundary
          disp_tmp((1+nn):(m+nn),j)=(min([nn+j-mm,2*nn+1-2*mm],[],2)*xDisp_mod(:,1,lln)+...
      	                             max([nn-j-mm,0],[],2)*xDisp_mod(:,n,lln))/(2*nn+1-2*mm);
          disp_tmp((1+nn):(m+nn),n+nn+j)=(max([j-mm,0],[],2)*xDisp_mod(:,1,lln)+...
	                                  min([2*nn+1-j-mm,2*nn+1-2*mm],[],2)*xDisp_mod(:,n,lln))/(2*nn+1-2*mm);
        end
        for j=1:1:nn
          disp_tmp(j,:)=(min([nn+j-mm,2*nn+1-2*mm],[],2)*disp_tmp(1+nn,:)+...
  	                 max([nn-j-mm,0],[],2)*disp_tmp(m+nn,:))/(2*nn+1-2*mm);
          disp_tmp(m+nn+j,:)=(max([j-mm,0],[],2)*disp_tmp(1+nn,:)+...
	                      min([2*nn+1-j-mm,2*nn+1-2*mm],[],2)*disp_tmp(m+nn,:))/(2*nn+1-2*mm);
        end
        disp_tmp2=fftshift(fft2(disp_tmp));
%        figure
%        subplot(1,2,1)
%        title(gcf,'real part');
%        plot3(xCoor_mod,yCoor_mod,real(disp_tmp2(1:m,1:n)))
%        subplot(1,2,2)
%        plot3(xCoor_mod,yCoor_mod,imag(disp_tmp2(1:m,1:n)))
%        title(gcf,'imaginary part');
        % get the values of the middle point and the margins for not zeroing the coefficients
        if mod(m,2)==0
          midm=m/2+1+nn;
        else
          midm=(m+1)/2+nn;
        end
        if mod(n,2)==0
          midn=n/2+1+nn;
        else
          midn=(n+1)/2+nn;
        end
        marm=min([nKeep,(midm-4)],[],2);
        marn=min([nKeep,(midn-4)],[],2);
        disp_tmp=zeros(m+2*nn,n+2*nn);
	disp_tmp(midm+[-marm-3:1:marm+3],midn+[-marn-3:1:marn+3])=0.25;
	disp_tmp(midm+[-marm-2:1:marm+2],midn+[-marn-2:1:marn+2])=0.5;
	disp_tmp(midm+[-marm-1:1:marm+1],midn+[-marn-1:1:marn+1])=0.75;
	disp_tmp(midm+[-marm:1:marm],midn+[-marn:1:marn])=1;
	disp_tmp=disp_tmp.*disp_tmp2;
        disp_tmp2=ifft2(ifftshift(disp_tmp));
        xDisp_mod(:,:,lln)=real(disp_tmp2((1+nn):(m+nn),(1+nn):(n+nn)));
      end
    elseif false% windowed fourier filter
      % separer le domain en sous-domaines aussi carres que possible,
      %  faire une analyse de fourier sur chaque sous-domaine,
      %  reconstruire sur chaque sous-domaine (cf fourier filer),
      %  et utiliser des fonctions chapeau pour reconstruire le domaine entier
    elseif false% wavelet filter
      % voir ce qui existe deja dans matlab
    end
    % check the noise on the x- and y-component of the displacement on the boundaries of the domain
    for lln=1:1:size(Ln,2)
      figure
      set(gcf,'color','w');
      subplot(1,2,1)
      hold on
      plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:),'k--');
      plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:)+xDisp_mod(1+bymin,:,lln),'b-');
      plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:),'k--');
      plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:)+xDisp_mod(floor(m/2),:,lln),'b-');
      plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:),'k--');
      plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:)+xDisp_mod(m-bymax,:,lln),'b-');
      plot(xCoor_mod(:,1+bxmin),yCoor_mod(:,1+bxmin),'k--');
      plot(xCoor_mod(:,1+bxmin)+xDisp_mod(:,1+bxmin,lln),yCoor_mod(:,1+bxmin),'r-');
      plot(xCoor_mod(:,floor(n/2)),yCoor_mod(:,floor(n/2)),'k--');
      plot(xCoor_mod(:,floor(n/2))+xDisp_mod(:,floor(n/2),lln),yCoor_mod(:,floor(n/2)),'r-');
      plot(xCoor_mod(:,n-bxmax),yCoor_mod(:,n-bxmax),'k--');
      plot(xCoor_mod(:,n-bxmax)+xDisp_mod(:,n-bxmax,lln),yCoor_mod(:,n-bxmax),'r-');
      title('Smoothed data: lateral disp. along some lines of the domain');
      legend('lines in the domain');
      axis equal
      subplot(1,2,2)
      hold on
      plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:),'k--');
      plot(xCoor_mod(1+bymin,:),yCoor_mod(1+bymin,:)+yDisp_mod(1+bymin,:,lln),'b-');
      plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:),'k--');
      plot(xCoor_mod(floor(m/2),:),yCoor_mod(floor(m/2),:)+yDisp_mod(floor(m/2),:,lln),'b-');
      plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:),'k--');
      plot(xCoor_mod(m-bymax,:),yCoor_mod(m-bymax,:)+yDisp_mod(m-bymax,:,lln),'b-');
      plot(xCoor_mod(:,1+bxmin),yCoor_mod(:,1+bxmin),'k--');
      plot(xCoor_mod(:,1+bxmin)+yDisp_mod(:,1+bxmin,lln),yCoor_mod(:,1+bxmin),'r-');
      plot(xCoor_mod(:,floor(n/2)),yCoor_mod(:,floor(n/2)),'k--');
      plot(xCoor_mod(:,floor(n/2))+yDisp_mod(:,floor(n/2),lln),yCoor_mod(:,floor(n/2)),'r-');
      plot(xCoor_mod(:,n-bxmax),yCoor_mod(:,n-bxmax),'k--');
      plot(xCoor_mod(:,n-bxmax)+yDisp_mod(:,n-bxmax,lln),yCoor_mod(:,n-bxmax),'r-');
      title('Axial disp. along some lines of the domain (if no disp., the dashed and colored lines superpose)');
      legend('Lines in the domain');
      axis equal
 %     figure
 %     subplot(1,2,1)
 %     surf(xCoor_mod,yCoor_mod,xDisp_mod(:,:,lln),'edgealpha',0);
 %     view(0,90);
 %     axis equal;
 %     title('Initial data: lateral displacement');
 %     subplot(1,2,2)
 %     surf(xCoor_mod,yCoor_mod,yDisp_mod(:,:,lln),'edgealpha',0);
 %     view(0,90);
 %     axis equal;
 %     title('axial displacement');
    end
    xDisp(ay2,ax2,Ln)=xDisp_mod;
    yDisp(ay2,ax2,Ln)=yDisp_mod;
  end% if false% use fitler
end

%if plotBCs  xDispNlms=xDisp; end

% check the incompressibility condition in the data given by tim. Maybe use a 3D median filter before differentiating... or stable numerical differentiation
%  see then where the incompressibility condition is worst violated; emphasis on the faces of the domain(as it is where the BCs will be imposed)
if false
% size of the data
ex = size(xDisp,2);
ey = size(xDisp,1);
ez = size(xDisp,3);
  disp('enforcing the incompressibily condition (assumes Uy is great quality)');
  
  % compute divu for the original displacements
  xDisp_smooth=create_input_medianFilter3D(xDisp(:,:,:,Ln));
  yDisp_smooth=yDisp(:,:,:,Ln);
  zDisp_smooth=create_input_medianFilter3D(zDisp(:,:,:,Ln));
  gradxXdisp=(xDisp_smooth(2:1:(ey-1),3:1:ex,2:1:(ez-1),:)-xDisp_smooth(2:1:(ey-1),1:1:(ex-2),2:1:(ez-1),:))/(xCoor(1,3,1)-xCoor(1,1,1));
  gradyYdisp=(yDisp_smooth(3:1:ey,2:1:(ex-1),2:1:(ez-1),:)-yDisp_smooth(1:1:(ey-2),2:1:(ex-1),2:1:(ez-1),:))/(yCoor(3,1,1)-yCoor(1,1,1));
  gradzZdisp=(zDisp_smooth(2:1:(ey-1),2:1:(ex-1),3:1:ez,:)-zDisp_smooth(2:1:(ey-1),2:1:(ex-1),1:1:(ez-2),:))/(zCoor(1,1,3)-zCoor(1,1,1));
  divu=gradxXdisp+gradyYdisp+gradzZdisp;
  % assumes Uz=cst and computes Ux accordingly.
  xDisp_theory=create_input_theoryPaul(yDisp(:,:,:,Ln),xDisp(:,:,:,Ln),xCoor(1,2,1)-xCoor(1,1,1),yCoor(2,1,1)-yCoor(1,1,1),zCoor(1,1,2)-zCoor(1,1,1));
  gradxXdisp_theory=(xDisp_theory(2:1:(ey-1),3:1:ex,2:1:(ez-1),:)-xDisp_theory(2:1:(ey-1),1:1:(ex-2),2:1:(ez-1),:))/(xCoor(1,3,1)-xCoor(1,1,1));
  divu2=gradxXdisp_theory+gradyYdisp;

  % iterative approach
  %[xDisp_theory3,zDisp_theory3]=create_input_theoryJFiterative(yDisp(:,:,:,Ln),xDisp(:,:,:,Ln),zDisp(:,:,:,Ln),xCoor(1,2,1)-xCoor(1,1,1),yCoor(2,1,1)-yCoor(1,1,1),zCoor(1,1,2)-zCoor(1,1,1));
  [xDisp_theory3,zDisp_theory3]=create_input_theoryJFiterative2(yDisp(:,:,:,Ln),xDisp(:,:,:,Ln),zDisp(:,:,:,Ln),xCoor(1,2,1)-xCoor(1,1,1),yCoor(2,1,1)-yCoor(1,1,1),zCoor(1,1,2)-zCoor(1,1,1));
  gradxXdisp_theory3=(xDisp_theory3(2:1:(ey-1),3:1:ex,2:1:(ez-1),:)-xDisp_theory3(2:1:(ey-1),1:1:(ex-2),2:1:(ez-1),:))/(xCoor(1,3,1)-xCoor(1,1,1));
  gradzZdisp_theory3=(zDisp_theory3(2:1:(ey-1),2:1:(ex-1),3:1:ez,:)-zDisp_theory3(2:1:(ey-1),2:1:(ex-1),1:1:(ez-2),:))/(zCoor(1,1,3)-zCoor(1,1,1));
  divu3=gradxXdisp_theory3+gradyYdisp+gradzZdisp_theory3;
  % iterative approach
  [xDisp_theory4,zDisp_theory4]=create_input_theory3(yDisp(:,:,:,Ln),xDisp_smooth,zDisp_smooth,xCoor(1,2,1)-xCoor(1,1,1),yCoor(2,1,1)-yCoor(1,1,1),zCoor(1,1,2)-zCoor(1,1,1));
%  [xDisp_theory4,zDisp_theory4]=create_input_theory2(yDisp(:,:,:,Ln),xDisp(:,:,:,Ln),zDisp(:,:,:,Ln),xCoor(1,2,1)-xCoor(1,1,1),yCoor(2,1,1)-yCoor(1,1,1),zCoor(1,1,2)-zCoor(1,1,1));
%  xDisp_theory4=create_input_medianFilter3D(xDisp_theory4);
%  zDisp_theory4=create_input_medianFilter3D(zDisp_theory4);
  gradxXdisp_theory4=(xDisp_theory4(2:1:(ey-1),3:1:ex,2:1:(ez-1),:)-xDisp_theory4(2:1:(ey-1),1:1:(ex-2),2:1:(ez-1),:))/(xCoor(1,3,1)-xCoor(1,1,1));
  gradzZdisp_theory4=(zDisp_theory4(2:1:(ey-1),2:1:(ex-1),3:1:ez,:)-zDisp_theory4(2:1:(ey-1),2:1:(ex-1),1:1:(ez-2),:))/(zCoor(1,1,3)-zCoor(1,1,1));
  divu4=gradxXdisp_theory4+gradyYdisp+gradzZdisp_theory4;
% natural colorscale
  [XX,YY,ZZ]=meshgrid(1:1:(ex-2),1:1:(ey-2),1:1:(ez-2));
  [XXa,YYa,ZZa]=meshgrid(1:1:ex,1:1:ey,1:1:ez);
 if false
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,divu,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('divergence u');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,divu,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]);
  title('divergence u (inside view)');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradxXdisp,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Ux,x');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradyYdisp,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Uy,y');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradzZdisp,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Uz,z');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
 end
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,divu2,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
%  title('divergence u theory form Paul');
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,divu2,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]);
%  title('divergence u theory from Paul');
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,gradxXdisp_theory,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
%  title('Ux,x theory from Paul');
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
if false
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,divu3,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('divergence u theory form JF');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,divu3,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]);
  title('divergence u theory from JF');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradxXdisp_theory3,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Ux,x theory from JF');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradzZdisp_theory3,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Uz,z theory from JF');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,xDisp_theory3(2:1:(ey-1),2:1:(ex-1),2:1:(ez-1)),[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Ux theory from JF');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,zDisp_theory3(2:1:(ey-1),2:1:(ex-1),2:1:(ez-1)),[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
  title('Uz theory from JF');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  end
  mar=1.0000001;
  figure
  set(gcf,'color','w');
  subplot(2,2,1); sl11=slice(XX,YY,ZZ,divu,[1+mar,ex-2-mar],[1+mar,ey-2-mar],[1+mar,ez-2-mar]); set(sl11,'edgealpha',0); title('div(u) original');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,2); sl12=slice(XX,YY,ZZ,divu2,[1+mar,ex-2-mar],[1+mar,ey-2-mar],[1+mar,ez-2-mar]); set(sl12,'edgealpha',0); title('div(u) Paul');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,3); sl13=slice(XX,YY,ZZ,divu3,[1+mar,ex-2-mar],[1+mar,ey-2-mar],[1+mar,ez-2-mar]); set(sl13,'edgealpha',0); title('div(u) JF iterative');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,4); sl14=slice(XX,YY,ZZ,divu4,[1+mar,ex-2-mar],[1+mar,ey-2-mar],[1+mar,ez-2-mar]); set(sl14,'edgealpha',0); title('div(u) theory');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');caxis([-0.05;0.05]);
  figure
  set(gcf,'color','w');
  subplot(2,2,1);sl21=slice(XX,YY,ZZ,divu,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]); set(sl21,'edgealpha',0); title('div(u) original');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,2);sl22=slice(XX,YY,ZZ,divu2,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]); set(sl22,'edgealpha',0); title('div(u) Paul');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,3);sl23=slice(XX,YY,ZZ,divu3,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]); set(sl23,'edgealpha',0); title('div(u) JF iterate');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,4);sl24=slice(XX,YY,ZZ,divu4,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]); set(sl24,'edgealpha',0); title('div(u) theory');
  title('divergence u (cut inside)');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');caxis([-0.05;0.05]);
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradxXdisp_theory4,[1+mar,ex-2-mar],[1+mar,ey-2-mar],[1+mar,ez-2-mar]);
  title('Ux,x theory');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,gradzZdisp_theory4,[1+mar,ex-2-mar],[1+mar,ey-2-mar],[1+mar,ez-2-mar]);
  title('Uz,z theory');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  subplot(2,2,1);sl51=slice(XXa,YYa,ZZa,xDisp_smooth,[1+mar,ex-mar],[1+mar,ey-mar],[1+mar,ez-mar]); set(sl51,'edgealpha',0); title('Ux original');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,2);sl52=slice(XXa,YYa,ZZa,xDisp_theory,[1+mar,ex-mar],[1+mar,ey-mar],[1+mar,ez-mar]); set(sl52,'edgealpha',0); title('Ux Paul');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,3);sl53=slice(XXa,YYa,ZZa,xDisp_theory3,[1+mar,ex-mar],[1+mar,ey-mar],[1+mar,ez-mar]); set(sl53,'edgealpha',0); title('Ux JF iterate');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,4);sl54=slice(XXa,YYa,ZZa,xDisp_theory4,[1+mar,ex-mar],[1+mar,ey-mar],[1+mar,ez-mar]); set(sl54,'edgealpha',0); title('Ux theory');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  toto=zeros(ey,ex,ez);
  subplot(2,2,1);sl61=slice(XXa,YYa,ZZa,zDisp_smooth,[1+mar,ex-mar],[1+mar,ey-mar],[1+mar,ez-mar]); set(sl61,'edgealpha',0); title('Uz original');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,2);sl62=slice(XXa,YYa,ZZa,toto,[1+mar,ex-0-mar],[1+mar,ey-0-mar],[1+mar,ez-0-mar]); set(sl62,'edgealpha',0); title('Uz Paul');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,3);sl63=slice(XXa,YYa,ZZa,zDisp_theory3,[1+mar,ex-0-mar],[1+mar,ey-0-mar],[1+mar,ez-0-mar]); set(sl63,'edgealpha',0); title('Uz JF iterate');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  subplot(2,2,4);sl64=slice(XXa,YYa,ZZa,zDisp_theory4,[1+mar,ex-0-mar],[1+mar,ey-0-mar],[1+mar,ez-0-mar]); set(sl64,'edgealpha',0); title('Uz theory');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
% different colorscale
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,divu,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
%  title('divergence u'); caxis([-0.05,0.05]);
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,divu,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]);
%  title('divergence u'); caxis([-0.05,0.05]);
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,divu2,[1.0000001,ex-2.0000001],[1.0000001,ey-2.0000001],[1.0000001,ez-2.0000001]);
%  title('divergence u theory form Paul'); caxis([-0.05,0.05]);
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
%  figure
%  set(gcf,'color','w');
%  slice(XX,YY,ZZ,divu2,[floor(ex/2)],[floor(ey/2)],[floor(ez/2)]);
%  title('divergence u theory from Paul'); caxis([-0.05,0.05]);
%  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');

  clear xDisp_smooth; clear yDisp_smooth; clear zDisp_smooth;
%  xDisp(:,:,:,Ln)=xDisp_theory;disp('modifying xDisp, assumes Uz=cst');
  %xDisp(:,:,:,Ln)=xDisp_theory3;disp('modifying xDisp');zDisp(:,:,:,Ln)=zDisp_theory3;disp('modifying zDisp');
  xDisp(:,:,:,Ln)=xDisp_theory4;disp('modifying xDisp');zDisp(:,:,:,Ln)=zDisp_theory4;disp('modifying zDisp');
  %return;
% apply the margins
%ex=ex-sideMargin;
ex=ex-lateMargin;
ey=ey-bottomMargin;
%ez=ez-sideMargin;
ez=ez-elevMargin;
end% if false

% Reducing xDisp and yDisp arrays
xDisp_mod=xDisp(ay,ax,az,Ln);
yDisp_mod=yDisp(ay,ax,az,Ln);
zDisp_mod=zDisp(ay,ax,az,Ln);
pDisp_mod=pDisp(ay,ax,az,Ln);
xCoor_mod=xCoor(ay,ax,az);
yCoor_mod=yCoor(ay,ax,az);
zCoor_mod=zCoor(ay,ax,az);
m=size(xDisp_mod,1);
n=size(xDisp_mod,2);
o=size(xDisp_mod,3);

% compute if 2 frames are different
if false
  tmp=yDisp_mod(:,:,:,1);
  normFrame1=norm(tmp(:),2);
  tmp=yDisp_mod(:,:,:,2);
  normFrame2=norm(tmp(:),2);
  normRef=normFrame2;
  for c=1.9:0.02:2.6
    tmp=c*yDisp_mod(:,:,:,1)-yDisp_mod(:,:,:,2);
    normtemp=norm(tmp(:),2);
    if normRef>normtemp
      normRef=normtemp;
    end
  end
  normRef^2/(normFrame1*normFrame2)
return
end% if false

% analysis of the data
if false
for lln=1:size(Ln,2)
  [XX,YY,ZZ]=meshgrid(1:n,1:m,1:o);% create a new grid
  % analysis of axial displacements
  tmp=yDisp_mod(1,:,:,lln);
  disp(['avg axial disp close to the transducer array: ' num2str(mean(tmp(:)))]);
  disp(['std axial disp close to the transducer array: ' num2str(std(tmp(:)))]);
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,yDisp_mod,[1.0000001,n-0.0000001],[1.0000001,m-0.0000001],[1.0000001,o-0.0000001]);
  title('axial displacements (axial direction = y-axis)');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  tmp=yDisp_mod(size(yDisp_mod,1),:,:,lln);
  disp(['avg axial disp far from the transducer array: ' num2str(mean(tmp(:)))]);
  disp(['std axial disp far from the transducer array: ' num2str(std(tmp(:)))]);

  % analysis of the lateral displacements
  tmp=xDisp_mod(:,1,:,lln);
  disp(['avg lateral disp on lateral face 1: ' num2str(mean(tmp(:)))]);
  disp(['std lateral disp on lateral face 1: ' num2str(std(tmp(:)))]);
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,xDisp_mod,[1.0000001,n-0.0000001],[1.0000001,m-0.0000001],[1.0000001,o-0.0000001]);
  title('lateral displacements (lateral direction = x-axis)');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,xDisp_mod,[1.0000001,n-0.0000001],[1.0000001,m-0.0000001],[1.0000001,o-0.0000001]);
  title('lateral displacements (lateral direction = x-axis) [new colorscale]');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)'); caxis([-0.8,0.2]);
  tmp=xDisp_mod(:,size(xDisp_mod,2),:,lln);
  disp(['avg lateral disp on lateral face 2: ' num2str(mean(tmp(:)))]);
  disp(['std lateral disp on lateral face 2: ' num2str(std(tmp(:)))]);

  % analysis of the elevational displacements
  tmp=zDisp_mod(:,:,1,lln);
  disp(['avg elevational disp on elevational face 1: ' num2str(mean(tmp(:)))]);
  disp(['std elevational disp on elevational face 1: ' num2str(std(tmp(:)))]);
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,zDisp_mod,[1.0000001,n-0.0000001],[1.0000001,m-0.0000001],[1.0000001,o-0.0000001]);
  title('elevational displacements (elevational direction = z-axis)');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)');
  figure
  set(gcf,'color','w');
  slice(XX,YY,ZZ,zDisp_mod,[1.0000001,n-0.0000001],[1.0000001,m-0.0000001],[1.0000001,o-0.0000001]);
  title('elevational displacements (elevational direction = z-axis) [new colorscale]');
  colorbar; xlabel('x (lateral)'); ylabel('y (axial)'); zlabel('z (elevational)'); caxis([0.1,0.7]);
  tmp=zDisp_mod(:,:,size(zDisp_mod,3),lln);
  disp(['avg elevational disp on elevational face 2: ' num2str(mean(tmp(:)))]);
  disp(['std elevational disp on elevational face 2: ' num2str(std(tmp(:)))]);
  clear XX; clear YY; clear ZZ;
end%for lln=1:size(Ln,2)
return
end% if false


% node_coor is an assembling of the node number and its x and y
% coordinates

max_node=m*n*o;
node_coor=zeros(max_node,4);
node_coor(:,1)=[1:1:max_node]';
nodesNum=zeros([n,m,o]);
nodesNum(:)=node_coor(:,1);

xCoor_mod2=zeros(n,m,o);% make it ndgrid format.
yCoor_mod2=zeros(n,m,o);
zCoor_mod2=zeros(n,m,o);
for k=1:1:o
  xCoor_mod2(:,:,k)=xCoor_mod(:,:,k)';
  yCoor_mod2(:,:,k)=yCoor_mod(:,:,k)';
  zCoor_mod2(:,:,k)=zCoor_mod(:,:,k)';
end
xCoor_mod=xCoor_mod2;
yCoor_mod=yCoor_mod2;
zCoor_mod=zCoor_mod2;
clear xCoor_mod2;clear yCoor_mod2;clear zCoor_mod2;
node_coor(:,2)=xCoor_mod(:);
node_coor(:,3)=yCoor_mod(:);
node_coor(:,4)=zCoor_mod(:);



% Element numbering and its respective nodes
if elType==32||elType==31||elType==305% making tetras (split every hexa in 5 tetras)
  elem_max=5*(n-1)*(m-1)*(o-1);
  elem=zeros(elem_max,5);
  elem(:,1)=[1:1:elem_max]';
  ind=[1:10:elem_max];% number of the elements built the same way
  tmp1=nodesNum(1:1:(n-1),1:1:(m-1),2:1:o);% node 1 of hexas ...
  tmp2=nodesNum(1:1:(n-1),2:1:m,2:1:o);
  tmp3=nodesNum(2:1:n,1:1:(m-1),2:1:o);
  tmp4=nodesNum(2:1:n,2:1:m,2:1:o);
  tmp5=nodesNum(1:1:(n-1),1:1:(m-1),1:1:(o-1));
  tmp6=nodesNum(1:1:(n-1),2:1:m,1:1:(o-1));
  tmp7=nodesNum(2:1:n,1:1:(m-1),1:1:(o-1));
  tmp8=nodesNum(2:1:n,2:1:m,1:1:(o-1));
  
  % the following formula would work only for even number of elements
  %oddIndex=[1:2:((n-1)*(m-1)*(o-1))]';
  %evenIndex=[2:2:((n-1)*(m-1)*(o-1))]';
  [numi,numj,numk]=ndgrid([1:1:(n-1)],[1:1:(m-1)],[1:1:(o-1)]);
  tmp=numi(:)+numj(:)+numk(:);
  oddIndex=find(mod(tmp,2)==1);
  evenIndex=find(mod(tmp,2)==0);

  elem(1+(oddIndex-1)*5,2:5)=[tmp1(oddIndex),tmp2(oddIndex),tmp3(oddIndex),tmp5(oddIndex)];
  elem(2+(oddIndex-1)*5,2:5)=[tmp2(oddIndex),tmp4(oddIndex),tmp3(oddIndex),tmp8(oddIndex)];
  elem(3+(oddIndex-1)*5,2:5)=[tmp2(oddIndex),tmp3(oddIndex),tmp5(oddIndex),tmp8(oddIndex)];
  elem(4+(oddIndex-1)*5,2:5)=[tmp5(oddIndex),tmp6(oddIndex),tmp2(oddIndex),tmp8(oddIndex)];
  elem(5+(oddIndex-1)*5,2:5)=[tmp5(oddIndex),tmp8(oddIndex),tmp3(oddIndex),tmp7(oddIndex)];
  elem(1+(evenIndex-1)*5,2:5)=[tmp1(evenIndex),tmp2(evenIndex),tmp4(evenIndex),tmp6(evenIndex)];
  elem(2+(evenIndex-1)*5,2:5)=[tmp1(evenIndex),tmp4(evenIndex),tmp3(evenIndex),tmp7(evenIndex)];
  elem(3+(evenIndex-1)*5,2:5)=[tmp1(evenIndex),tmp6(evenIndex),tmp4(evenIndex),tmp7(evenIndex)];
  elem(4+(evenIndex-1)*5,2:5)=[tmp1(evenIndex),tmp6(evenIndex),tmp7(evenIndex),tmp5(evenIndex)];
  elem(5+(evenIndex-1)*5,2:5)=[tmp6(evenIndex),tmp4(evenIndex),tmp7(evenIndex),tmp8(evenIndex)];
elseif elType==105% make hexas
  elem_max=(n-1)*(m-1)*(o-1);
  elem=zeros(elem_max,9);
  elem(:,1)=[1:1:elem_max]';
  tmp=nodesNum(1:1:(n-1),1:1:(m-1),2:1:o);elem(:,2)=tmp(:);
  tmp=nodesNum(1:1:(n-1),2:1:m,2:1:o);elem(:,3)=tmp(:);
  tmp=nodesNum(2:1:n,1:1:(m-1),2:1:o);elem(:,4)=tmp(:);
  tmp=nodesNum(2:1:n,2:1:m,2:1:o);elem(:,5)=tmp(:);
  tmp=nodesNum(1:1:(n-1),1:1:(m-1),1:1:(o-1));elem(:,6)=tmp(:);
  tmp=nodesNum(1:1:(n-1),2:1:m,1:1:(o-1));elem(:,7)=tmp(:);
  tmp=nodesNum(2:1:n,1:1:(m-1),1:1:(o-1));elem(:,8)=tmp(:);
  tmp=nodesNum(2:1:n,2:1:m,1:1:(o-1));elem(:,9)=tmp(:);
else
  disp('elType has a value that is not recognized ... exiting');
  return;
end

% check if the displacements will lead to an inverted mesh
if false
end% if false

% check the volume of a series of frames
if false
  volumes=zeros(size(Ln,2)+1,1);% total volume
  elemVols=zeros(size(Ln,2)+1,size(elem,1));% volumes of each element
  for i=1:1:size(elem,1)
    v1=node_coor(elem(i,3),2:4)-node_coor(elem(i,2),2:4);
    v2=node_coor(elem(i,4),2:4)-node_coor(elem(i,2),2:4);
    v3=node_coor(elem(i,5),2:4)-node_coor(elem(i,2),2:4);
    volu=dot(cross(v2,v1),v3)/6;
    elemVols(1,i)=volu;
    volumes(1,1)=volumes(1,1)+volu;
  end% for i=1:1:size(elem,1)
  disp(['the volume of the reference frame (diplacements=0) is: ' num2str(volumes(1,1))]);
  for lln=1:1:size(Ln,2)
    xCoor_mod2=zeros(n,m,o);% make it ndgrid format.% temporary variable
    yCoor_mod2=zeros(n,m,o);
    zCoor_mod2=zeros(n,m,o);
    for k=1:1:o
      xCoor_mod2(:,:,k)=xDisp_mod(:,:,k,lln)';
      yCoor_mod2(:,:,k)=yDisp_mod(:,:,k,lln)';
      zCoor_mod2(:,:,k)=zDisp_mod(:,:,k,lln)';
    end
    xCoor_mod=xCoor_mod2;
    yCoor_mod=yCoor_mod2;
    zCoor_mod=zCoor_mod2;
    node_coor2=zeros(size(xCoor_mod(:),1),3);
    node_coor2(:,2)=xCoor_mod(:)+node_coor(:,2);
    node_coor2(:,3)=yCoor_mod(:)+node_coor(:,3);
    node_coor2(:,4)=zCoor_mod(:)+node_coor(:,4);
    nNegVol=0;
    for i=1:1:size(elem,1)
      v1=node_coor2(elem(i,3),2:4)-node_coor2(elem(i,2),2:4);
      v2=node_coor2(elem(i,4),2:4)-node_coor2(elem(i,2),2:4);
      v3=node_coor2(elem(i,5),2:4)-node_coor2(elem(i,2),2:4);
      volu=dot(cross(v2,v1),v3)/6;
      if volu<0
        %disp(['warning: volume of element ' num2str(i) ' is negative']);
	nNegVol=nNegVol+1;
      end
      elemVols(1+lln,i)=volu;
      volumes(1+lln,1)=volumes(1+lln,1)+volu;
    end% for i=1:1:size(elem,1)
    disp(['the volume of frame ' num2str(Ln(lln)) ' is: ' num2str(volumes(1+lln,1))]);
    disp(['  ' num2str(nNegVol) ' elements have a negative volume in the deformed state']);
    %compute the histogram of the volume changes for all the elements
    elemVolChange=elemVols(1+lln,:)-elemVols(1,:);
    nBinsHisto=7;
    binBounds=[-1.0e-2,-4.0e-3,-1.0e-3,1.0e-3,4.0e-3,1.0e-2];
    histo=zeros(nBinsHisto,1);
    histo(1,1)=size(find(elemVolChange<binBounds(1,1)),2);
    for i=2:1:(nBinsHisto-1)
      histo(i,1)=size(find((elemVolChange<binBounds(1,i)).*(elemVolChange>binBounds(1,i-1))),2);
    end
    histo(nBinsHisto,1)=size(find(elemVolChange>binBounds(1,nBinsHisto-1)),2);
    disp('  histogram of the elemental volume change:');
    disp(['    ' num2str(histo')]);
    disp(['    (bin bounds: ' num2str(binBounds) ')']);
  end% for lln=1:size(Ln,2)
  return
end% if false

% enforce a constant volume by condencing the lateral displacements around their average value and also decrease the variations of the elevational displacements
if false
  volumes=zeros(size(Ln,2)+1,1);
  for i=1:1:size(elem,1)
    v1=node_coor(elem(i,3),2:4)-node_coor(elem(i,2),2:4);
    v2=node_coor(elem(i,4),2:4)-node_coor(elem(i,2),2:4);
    v3=node_coor(elem(i,5),2:4)-node_coor(elem(i,2),2:4);
    volu=dot(cross(v2,v1),v3)/6;
    volumes(1,1)=volumes(1,1)+volu;
  end% for i=1:1:size(elem,1)
  disp(['the volume of the reference frame (diplacements=0) is: ' num2str(volumes(1,1))]);
  for lln=1:1:size(Ln,2)
    xCoor_mod2=zeros(n,m,o);% make it ndgrid format.% temporary variable
    yCoor_mod2=zeros(n,m,o);
    zCoor_mod2=zeros(n,m,o);
    for k=1:1:o
      xCoor_mod2(:,:,k)=xDisp_mod(:,:,k,lln)';
      yCoor_mod2(:,:,k)=yDisp_mod(:,:,k,lln)';
      zCoor_mod2(:,:,k)=zDisp_mod(:,:,k,lln)';
    end
    node_coor2(:,2)=xCoor_mod2(:)+node_coor(:,2);
    node_coor2(:,3)=yCoor_mod2(:)+node_coor(:,3);
    node_coor2(:,4)=zCoor_mod2(:)+node_coor(:,4);
    for i=1:1:size(elem,1)
      v1=node_coor2(elem(i,3),2:4)-node_coor2(elem(i,2),2:4);
      v2=node_coor2(elem(i,4),2:4)-node_coor2(elem(i,2),2:4);
      v3=node_coor2(elem(i,5),2:4)-node_coor2(elem(i,2),2:4);
      volu=dot(cross(v2,v1),v3)/6;
      if volu<0
        disp(['warning: volume of element ' num2str(i) ' is negative']);
      end
      volumes(1+lln,1)=volumes(1+lln,1)+volu;
    end% for i=1:1:size(elem,1)
    disp(['the initial volume of frame ' num2str(Ln(lln)) ' is: ' num2str(volumes(1+lln,1))]);
    iter=0;
    meanLat=mean(xCoor_mod2(:));
    meanEle=mean(zCoor_mod2(:));
    volmem=volumes(1+lln,1);
    while ((abs(volumes(1+lln,1)-volumes(1,1))/volumes(1,1)>0.0004)&&(iter<50))||...
          ((volmem-volumes(1,1))*(volumes(1+lln,1)-volumes(1,1))<0)
      iter=iter+1;
      scaleFactor=1+((volumes(1,1)/volumes(1+lln,1))-1)*5;
      if ((volmem-volumes(1,1))*(volumes(1+lln,1)-volumes(1,1))<0)
        if volmem-volumes(1,1)<0
	  if scaleFactor>0.97
	    scaleFactor=0.97;
	  end
	else
	  if scaleFactor<1.03
	    scaleFactor=1.03;
	  end
	end
      end
      if scaleFactor>1.3
        scaleFactor=1.3;
      end
      if scaleFactor<0.75;
        scaleFactor=0.75;
      end
      disp(['scaleFactor = ' num2str(scaleFactor)]);
      %meanLat
      %disp(['lat. disp: ' num2str(min(xCoor_mod2(:))) ' ' num2str(max(xCoor_mod2(:)))]);
      xCoor_mod2=meanLat+(xCoor_mod2-meanLat)*scaleFactor;
      %disp(['lat. disp: ' num2str(min(xCoor_mod2(:))) ' ' num2str(max(xCoor_mod2(:)))]);
      %meanEle
      %disp(['elev. disp: ' num2str(min(zCoor_mod2(:))) ' ' num2str(max(zCoor_mod2(:)))]);
      zCoor_mod2=meanEle+(zCoor_mod2-meanEle)*(1+(scaleFactor-1)/2);
      %disp(['elev. disp: ' num2str(min(zCoor_mod2(:))) ' ' num2str(max(zCoor_mod2(:)))]);
      node_coor2(:,2)=xCoor_mod2(:)+node_coor(:,2);
      node_coor2(:,3)=yCoor_mod2(:)+node_coor(:,3);
      node_coor2(:,4)=zCoor_mod2(:)+node_coor(:,4);
      volmem2=volumes(1+lln,1);
      volumes(1+lln,1)=0;
      for i=1:1:size(elem,1)
        v1=node_coor2(elem(i,3),2:4)-node_coor2(elem(i,2),2:4);
        v2=node_coor2(elem(i,4),2:4)-node_coor2(elem(i,2),2:4);
        v3=node_coor2(elem(i,5),2:4)-node_coor2(elem(i,2),2:4);
        volu=dot(cross(v2,v1),v3)/6;
        %if volu<0
        %  disp(['warning: volume of element ' num2str(i) ' is negative']);
        %end
        volumes(1+lln,1)=volumes(1+lln,1)+volu;
      end% for i=1:1:size(elem,1)
      disp(['iter = ' num2str(iter) ', volume = ' num2str(volumes(1+lln,1)) ', jump = ' num2str(volumes(1+lln,1)-volmem2)]);
    end
    % get the displacements back into xDisp_mod, yDisp_mod and zDisp_mod
    xCoor_mod3=zeros(m,n,o);% make it ndgrid format.
    zCoor_mod3=zeros(m,n,o);
    for k=1:1:o
      xCoor_mod3(:,:,k)=xCoor_mod2(:,:,k)';
      zCoor_mod3(:,:,k)=zCoor_mod2(:,:,k)';
    end
    xDisp_mod(:,:,:,lln)=xCoor_mod3;
    zDisp_mod(:,:,:,lln)=zCoor_mod3;
    
    disp(['the final volume of frame ' num2str(Ln(lln)) ' is: ' num2str(volumes(1+lln,1)) '(iter=' num2str(iter) ')']);
  end% for lln=1:size(Ln,2)
end% if false

if addElem307
  disp(['[xmin,xmax,ymin,ymax,zmin,zmax]=[' num2str(min(xCoor_mod(:))) ',' num2str(max(xCoor_mod(:))) ',' num2str(min(yCoor_mod(:))) ',' num2str(max(yCoor_mod(:))) ',' num2str(min(zCoor_mod(:))) ',' num2str(max(zCoor_mod(:))) ']']);
  %nodesWithForce=nodesNum(1,:,:);% X=Xmin
  %nodesWithForce=nodesNum(n,:,:);% X=Xmax
  %nodesWithForce=nodesNum(:,1,:);disp(['chosen value of y=ymax=' num2str(yCoor_mod(1,1,1))]);% Y=Ymax
  nodesWithForce=nodesNum(:,m,:);disp(['chosen value of y=ymin=' num2str(yCoor_mod(1,m,1))]);% Y=Ymin
  indElemsWithForce=find(sum(ismember(elem(:,2:5),nodesWithForce(:)),2)>=0.5);% all the elements that have a node on the surface where the force is computed
  nElem307=size(indElemsWithForce,1);
  surfaceActive=(max(xCoor_mod(:))-min(xCoor_mod(:)))*(max(zCoor_mod(:))-min(zCoor_mod(:)));
  disp(['surfaceActive = ' num2str(surfaceActive)]);
end

if addElem611
  if notop
    nElem611 = 2 * 2 * ( (m-1)*(n-1) + (m-1)*(o-1) + (n-1)*(o-1) ) ...
        - 2*((n-1)*(o-1)) - 4*((n-1) + (o-1)) ; 
  elseif notopandbottom
    nElem611 = 2 * 2 * ( (m-1)*(n-1) + (m-1)*(o-1) + (n-1)*(o-1) ) ...
        - 4*((n-1)*(o-1)) - 8*((n-1) + (o-1)) ; 
  else
    nElem611 = 2 * 2 * ( (m-1)*(n-1) + (m-1)*(o-1) + (n-1)*(o-1) ) ;
  end
else
  nElem611=0;
end

% Creates a file called input.m (lots of comments are printed for clarity)
if addElem611
  fid = fopen(['with611-' filename(1:(end-4)) '.in'], 'wt');
else
  fid = fopen([filename(1:(end-4)) '.in'], 'wt');
end
fprintf(fid,['%% this is an input file for NLACE created on ' datestr(clock) '\n']);
fprintf(fid,['%%  from ' filename '\n']);
fprintf(fid,'\n');
fprintf(fid,'%% NLACE is based on an optimization algorithm (L-BFGS-B, ...)\n');
fprintf(fid,'%%  and a solver (Non-Linear elastic solver). NLACE can solve\n');
fprintf(fid,'%%  linear inverse problems and forward problems as well.\n');
fprintf(fid,'\n');
if size(Ln,2)==1
  fprintf(fid,['%% In this input file, there is ' num2str(size(Ln,2)) ' measurement of displacements\n']);
  fprintf(fid,['%% It is compression frame #' num2str(Ln(1)) '\n']);
else
  fprintf(fid,['%% In this input file, there are ' num2str(size(Ln,2)) ' different measurements of displacements\n']);
  fprintf(fid,['%% They are the compression frames #' num2str(Ln) '\n']);
end
fprintf(fid,['%%  with options: reduce=' num2str(reduce) ', topMargin=' num2str(topMargin) '\n']);
%fprintf(fid,['%%  bottomMargin=' num2str(bottomMargin) ', sideMargin=' num2str(sideMargin) '\n']);
fprintf(fid,['%%  bottomMargin=' num2str(bottomMargin) ', lateMargin=' num2str(lateMargin) ', elevMargin=' num2str(elevMargin) '\n']);
fprintf(fid,'\n');
fprintf(fid,'%% Every section begins with a keyword and the order of the sections/keywords does not\n');
fprintf(fid,'%%  matter. The keywords are words beginning with: ''femo'', ''opto'', ''solo'', ''outp'',\n');
fprintf(fid,'%%  ''opts'', ''timi'', ''pbty'', ''optp'', ''end '', ''elem'', ''optb'', ''datn'', ''coor'',\n');
fprintf(fid,'%%  ''conn'', ''mate'', ''boun'', ''rest'', ''date'' or ''meas''.\n'); 
fprintf(fid,'\n');
fprintf(fid,'%% Comments can be done between any section but no comment is allowed between the keyword\n');
fprintf(fid,'%%  and its attached data. It is recommended to use ''%%'' for comments but it is\n');
fprintf(fid,'%%  not mandatory as any line beginning with letters different from the above-mentionned\n');
fprintf(fid,'%%  keywords will be assumed to be a comment.\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% TYPE OF PROBLEM (OPTIONAL)\n');
fprintf(fid,'%% - pbtype: 1-> inverse problem (default); 2-> direct solve\n');
fprintf(fid,'pbtype\n');
if elType==305
  fprintf(fid,'%d\n',1);
elseif elType==31||elType==32
  fprintf(fid,'%d\n',2);
else
  disp('weird elType ... exiting');
  return
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% PRINTING TIMING INFO (FOR PROFILING - OPTIONAL)\n');
fprintf(fid,'%%timingPrint\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
if elType==305
fprintf(fid,'%% OPTIONS FOR THE BOUND-CONSTRAINED OPTIMIZATION (or BOX-CONSTRAINED)\n');
fprintf(fid,'%%  - iopt   : optimization strategy (2->L-BFGS-B; 4->ASA_CG; 5->gencan)\n');
fprintf(fid,'%%  - niter  : maximum number of iterations\n');
fprintf(fid,'%%  - bfgsM  : value of parameter m in L-BFGS-B (usually 10<m<25)\n');
fprintf(fid,'%%  - noutput: interval for saving data\n');
fprintf(fid,'%%  - mnty   : maximum number of fields (types of variables) to be optimized (cf optbounds)\n');
fprintf(fid,'%% iopt  niter  bfgsM  noutput  mnty\n');
fprintf(fid,'optoptions\n');
fprintf(fid,'%d %d %d %d %d\n',2,2000,10,500,2);
fprintf(fid,'\n');
fprintf(fid,'\n');
end
fprintf(fid,'%% OPTIONS FOR THE FEM CODE\n');
fprintf(fid,'%%  - nelem  : number of elements for the solver\n');
fprintf(fid,'%%  - npoin  : number of points for the solver\n');
fprintf(fid,'%%  - ndime  : dimension of the space (2D, 3D, ...)\n');
fprintf(fid,'%%  - mnode  : maximum number of nodes per element\n');
fprintf(fid,'%%  - mdofn  : maximum number of DOFs per node\n');
fprintf(fid,'%%  - nmat   : number of materials\n');
fprintf(fid,'%%  - mprops : maximum number of properties per material\n');
fprintf(fid,'%%  - nmeas  : number of measured displacement cases\n');
fprintf(fid,'%%  - mpoinbc: number of Dirichlet points in BC\n');
fprintf(fid,'%%  - tol    : tolerance for the convergence of the non-linear solver\n');
fprintf(fid,'%%  - lsteps : number of steps for loading the material (any value is ok)\n');
fprintf(fid,'%%  - ncontin: number of steps for transitioning material properties (any value is ok)\n');
fprintf(fid,'%% nelem  npoin  ndime  mnode  mdofn  nmat  mprops  nmeas  mpoinbc  tol  lsteps  ncontin\n');
fprintf(fid,'femoptions\n');
% modify to handle multiple frames
if elType==305
  if addElem307
    info1 = [5*(m-1)*(n-1)*(o-1)+nElem307 max_node 3 4 4 1+nElem307 10 size(Ln,2) 2*m*n+2*m*(o-2)+2*(n-2)*(o-2) 1.0e-11 5 1];
    
  elseif addElem611
    info1 = [5*(m-1)*(n-1)*(o-1)+nElem611 max_node 3 4 4 2 10 size(Ln,2) 2*m*n+2*m*(o-2)+2*(n-2)*(o-2) 1.0e-11 5 1];
  else
    info1 = [5*(m-1)*(n-1)*(o-1) max_node 3 4 4 1 10 size(Ln,2) 2*m*n+2*m*(o-2)+2*(n-2)*(o-2) 1.0e-11 5 1];
  end
elseif elType==31||elType==32
  info1 = [5*(m-1)*(n-1)*(o-1) max_node 3 4 3 1 5 1 0 1.0e-8 1 1];
else%if elType==105
  %info1 = [    (m-1)*(n-1) max_node   2     4     2     1     5  size(Ln,2)  2*(m+n)-4 1.0e-11  5   1];
  return;
end
fprintf(fid,'%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% OPTIONS FOR THE SOLVER (OPTIONAL)\n');
fprintf(fid,'%%  - nCores  : number of cores used by the solver (default nCores=1)\n');
fprintf(fid,'%%  - solveMethod: 1-> direct solver (default); 2-> iterative solver\n');
fprintf(fid,'%% nCores  solveMethod\n');
fprintf(fid,'soloptions\n');
fprintf(fid,'%d   %d\n',8,1);
fprintf(fid,'\n');
fprintf(fid,'\n');
if elType==305
fprintf(fid,'%% BOUNDS FOR THE BOUND-CONSTRAINED OPTIMIZATION (or BOX-CONSTRAINED)\n');
fprintf(fid,'%%  - lpa1: lower bound for parameter 1\n');
fprintf(fid,'%%  - upa1: upper bound for parameter 1\n');
fprintf(fid,'%%  - lpa2: lower bound for parameter 2\n');
fprintf(fid,'%%  - upa2: upper bound for parameter 2, etc\n');
fprintf(fid,'%% lpa1  upa1\n');
fprintf(fid,'%% lpa2  upa2 (there should be mnty lines)\n');
fprintf(fid,'optbounds\n');
fprintf(fid,'%f %f\n',0.01,4.0);
fprintf(fid,'%f %f\n',1.0,20.0);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% RESTART OPTION (commented by default)\n');
fprintf(fid,'%%restart\n');
fprintf(fid,'uncommentTheRestartKeywordAndPutYourFileNameHere\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
end
fprintf(fid,'%% NAME OF THE OUTPUT FILES (PREFIX ONLY - OPTIONAL)\n');
fprintf(fid,'output\n');
if addElem611
  fprintf(fid,['with611-' filename(1:(end-4)) '.out\n']);
else
  fprintf(fid,[filename(1:(end-4)) '.out\n']);
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% MATERIAL PROPERTIES\n');
if elType==305
  fprintf(fid,'%%  nGausPtsPerDirection regularizationType nonLinParamReg shearModReg cstTVD cstTVD2 T1 T2 stabilizationFactor activeStabilization\n');
  if addElem611
    fprintf(fid,'%%  nGaussPtsPerDir xk11 xk12 xk13 xk21 xk22 xk23 xk31 xk32 xk33\n');
  end
  fprintf(fid,'materialProperties\n');
  info3=[4  21  0.00000 0.00001  7.0d-03 7.0e-3 1.0 1.0 1.0 1];
  fprintf(fid,' %d %d %e %e %e %e %e %e %e %d\n',info3);
  if addElem611
    info32=[3 1.0d+02 0.000 0.000 0.000 5.0d+04 0.000 0.000 0.000 1.0d+02];
    fprintf(fid,' %d %e %e %e %e %e %e %e %e %e\n' ,info32);
  end
elseif elType==31||elType==32
  fprintf(fid,'%% nGausPtsPerDirection  alpha  dd(1)  dd(2)  dd(3)\n');
  fprintf(fid,'materialProperties\n');
  info3=[4   100 1 10 1];
  fprintf(fid,'%d %e %e %e %e\n',info3);
else%if elType==606
  %fprintf(fid,'%%  nGaussPtsPerDirection typeRegularisation gammaReg shearModReg cstTVD\n');
  %fprintf(fid,'materialProperties\n');
  %info3=[3  2  0.000 0.00001 1.0e-004];
  %fprintf(fid,' %d %d %f %f %f\n',info3);
  return;
end
if addElem307
  for i=1:1:nElem307
    if ismember(elem(indElemsWithForce(i),2),nodesWithForce)
      e1=elem(indElemsWithForce(i),2);
    else
      e1=0;
    end
    if ismember(elem(indElemsWithForce(i),3),nodesWithForce)
      e2=elem(indElemsWithForce(i),3);
    else
      e2=0;
    end
    if ismember(elem(indElemsWithForce(i),4),nodesWithForce)
      e3=elem(indElemsWithForce(i),4);
    else
      e3=0;
    end
    if ismember(elem(indElemsWithForce(i),5),nodesWithForce)
      e4=elem(indElemsWithForce(i),5);
    else
      e4=0;
    end
    info4=[4,e1,e2,e3,e4,totalForce*surfaceActive/surfaceForce,1.0];
    fprintf(fid,' %d %d %d %d %d %f %f\n',info4);
  end
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% ELEMENT SETS\n');
fprintf(fid,'%% initialElement  finalElement  elementType  numberOfNodesPerElement\n');
fprintf(fid,'elementSet(s)\n');
if elType==305
  info2=[1 5*(m-1)*(n-1)*(o-1) 305 4];
  fprintf(fid, '%d %d %d %d\n',info2);
elseif elType==31
  info2=[1 5*(m-1)*(n-1)*(o-1) 31 4];
  fprintf(fid, '%d %d %d %d\n',info2);
elseif elType==32
  info2=[1 5*(m-1)*(n-1)*(o-1) 32 4];
  fprintf(fid, '%d %d %d %d\n',info2);
else%if elType==606
  %info2=[1 (m-1)*(n-1) 606 4];
  %fprintf(fid, '%d %d %d %d\n',info2);
  return
end
if addElem307
  for i=1:1:nElem307
    info5=[5*(m-1)*(n-1)*(o-1)+i,5*(m-1)*(n-1)*(o-1)+i,307,4];
    fprintf(fid,'%d %d %d %d\n',info5);
  end
end
if addElem611
  info5=[ (5*(m-1)*(n-1)*(o-1) + 1) (5*(m-1)*(n-1)*(o-1) + nElem611) 611 3];
  fprintf(fid, '%d %d %d %d\n',info5);
end


fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% NODES COORDINATES\n');
fprintf(fid,'%% nodeNumber  x-coordinate  y-coordinate  z-coordinate\n');
fprintf(fid,'coordinates\n');
fprintf(fid,'%d %12.6e %12.6e %12.6e\n',node_coor');
fprintf(fid,'\n');
fprintf(fid,'\n');
% Writing connectivity data
fprintf(fid,'%% ELEMENTS TOPOLOGY\n');
fprintf(fid,'%% elementNumber nodeNumbers ...\n');
fprintf(fid,'connectivity\n');
if size(elem,2)==4
  fprintf(fid,' %d   %d   %d   %d\n',elem');
elseif size(elem,2)==5
  fprintf(fid,' %d   %d   %d   %d   %d\n',elem');
else
  disp('create_input.m, pb ... exiting');
  return;
end
if addElem307
  for i=1:1:nElem307
    elem307=[size(elem,1)+i,elem(indElemsWithForce(i),2:5)];
    fprintf(fid,'  %d   %d   %d   %d   %d\n',elem307);
  end
end
if addElem611
  %initialize some variables
  sheet = m*n;
  enum = 1;
  econn = zeros(nElem611,4);
    
  %top
  if notop==false&&notopandbottom==false 
    mtop = zeros(n,o); 
    iflip = false ;
    rflip = false ;
    for icol=1:o
      mtop(:, icol) = (icol-1)*sheet + sort([1:n],2,'descend')';
    end
    for irow = 1:(n-1)
      for iblock=1:(o-1)
        esquare = zeros(2,2);
        esquare(1,1) = mtop(irow,iblock);
        esquare(1,2) = mtop(irow,iblock + 1);
        esquare(2,1) = mtop(irow + 1, iblock);
        esquare(2,2) = mtop(irow + 1, iblock + 1);
        if iflip
          econn(enum,1) = enum;
          %element 1
          econn(enum,2) = esquare(2,2);
          econn(enum,3) = esquare(1,1);
          econn(enum,4) = esquare(2,1);
          enum = enum + 1;
          %element 2
          econn(enum,1) = enum;
          econn(enum,2) = esquare(1,1);
          econn(enum,3) = esquare(2,2);
          econn(enum,4) = esquare(1,2);
          enum = enum + 1;
          %flip the block
          iflip = not(iflip);
        else
          econn(enum,1) = enum;
          %element 1
          econn(enum,2) = esquare(2,1);
          econn(enum,3) = esquare(1,2);
          econn(enum,4) = esquare(1,1);
          enum = enum + 1;
          %element 2
          econn(enum,1) = enum;
          econn(enum,2) = esquare(1,2);
          econn(enum,3) = esquare(2,1);
          econn(enum,4) = esquare(2,2);
          enum = enum + 1;
          %flip the block
          iflip = not(iflip);
        end
      end
      rflip = not(rflip);
      iflip = rflip;
    end
  %else 
  end
  if notop == false 
    start = 1 ;
  else 
    start = 2;
  end   
  if notopandbottom == false
    dr = 0;
  else 
    dr = 1;
    start = 2;
  end
  %bottom
  if notopandbottom == false 
    mbottom = zeros(n,o); 
    iflip = false ;
    rflip = false ;
    for icol = 1:o
      mbottom(:,icol) =  (o-(icol-1))*sheet - (n-1) + sort([0:(n-1)],2,'descend') ; 
    end
    for irow = 1:(n-1)
      for iblock=1:(o-1)
        esquare = zeros(2,2);
        esquare(1,1) = mbottom(irow,iblock);
        esquare(1,2) = mbottom(irow,iblock + 1);
        esquare(2,1) = mbottom(irow + 1, iblock);
        esquare(2,2) = mbottom(irow + 1, iblock + 1);
        if iflip
          econn(enum,1) = enum;
          %element 1
          econn(enum,2) = esquare(2,2);
          econn(enum,3) = esquare(1,1);
          econn(enum,4) = esquare(2,1);
          enum = enum + 1;
          %element 2
          econn(enum,1) = enum;
          econn(enum,2) = esquare(1,1);
          econn(enum,3) = esquare(2,2);
          econn(enum,4) = esquare(1,2);
          enum = enum + 1;
          %flip the block
          iflip = not(iflip);
        else
          econn(enum,1) = enum;
          %element 1
          econn(enum,2) = esquare(2,1);
          econn(enum,3) = esquare(1,2);
          econn(enum,4) = esquare(1,1);
          enum = enum + 1;
          %element 2
          econn(enum,1) = enum;
          econn(enum,2) = esquare(1,2);
          econn(enum,3) = esquare(2,1);
          econn(enum,4) = esquare(2,2);
          enum = enum + 1;
          %flip the block
          iflip = not(iflip);
        end
      end
      rflip = not(rflip);
      iflip = rflip;
    end
  %else
  end
  %front
  mfront = zeros(m,o); 
  iflip = false ;
  rflip = false ;
  for icol = 1:o
    mfront(:,o - (icol-1)) = (icol - 1 )*sheet + [n:n:m*n]';
  end
  for irow = start:(m-1-dr)
    for iblock=1:(o-1)
      esquare = zeros(2,2);
      esquare(1,1) = mfront(irow,iblock);
      esquare(1,2) = mfront(irow,iblock + 1);
      esquare(2,1) = mfront(irow + 1, iblock);
      esquare(2,2) = mfront(irow + 1, iblock + 1);
      if iflip
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,2);
        econn(enum,3) = esquare(1,1);
        econn(enum,4) = esquare(2,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,1);
        econn(enum,3) = esquare(2,2);
        econn(enum,4) = esquare(1,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      else
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,1);
        econn(enum,3) = esquare(1,2);
        econn(enum,4) = esquare(1,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,2);
        econn(enum,3) = esquare(2,1);
        econn(enum,4) = esquare(2,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      end
    end
    rflip = not(rflip);
    iflip = rflip;
  end
  %back
  mback = zeros(m,o); 
  iflip = false ;
  rflip = false ;
  for icol = 1:o
    mback(:,icol) = (icol - 1 )*sheet + 1 + [0:n:n*n]';
  end
  for irow = start:(m-1-dr)
    for iblock=1:(o-1)
      esquare = zeros(2,2);
      esquare(1,1) = mback(irow,iblock);
      esquare(1,2) = mback(irow,iblock + 1);
      esquare(2,1) = mback(irow + 1, iblock);
      esquare(2,2) = mback(irow + 1, iblock + 1);
      if iflip
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,2);
        econn(enum,3) = esquare(1,1);
        econn(enum,4) = esquare(2,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,1);
        econn(enum,3) = esquare(2,2);
        econn(enum,4) = esquare(1,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      else
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,1);
        econn(enum,3) = esquare(1,2);
        econn(enum,4) = esquare(1,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,2);
        econn(enum,3) = esquare(2,1);
        econn(enum,4) = esquare(2,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      end
    end
    rflip = not(rflip);
    iflip = rflip;
  end
  %right
  mright = zeros(m,n); 
  iflip = false ;
  rflip = false ;
  for icol = 1:n
    mright(:,n - (icol - 1)) = (icol*ones(m,1) + [0:n:n*n]') ;
  end
  for irow = start:(m-1-dr)
    for iblock=1:(n-1)
      esquare = zeros(2,2);
      esquare(1,1) = mright(irow,iblock);
      esquare(1,2) = mright(irow,iblock + 1);
      esquare(2,1) = mright(irow + 1, iblock);
      esquare(2,2) = mright(irow + 1, iblock + 1);
      if iflip
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,2);
        econn(enum,3) = esquare(1,1);
        econn(enum,4) = esquare(2,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,1);
        econn(enum,3) = esquare(2,2);
        econn(enum,4) = esquare(1,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      else
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,1);
        econn(enum,3) = esquare(1,2);
        econn(enum,4) = esquare(1,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,2);
        econn(enum,3) = esquare(2,1);
        econn(enum,4) = esquare(2,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      end
    end
    rflip = not(rflip);
    iflip = rflip;
  end
  %left
  mleft = zeros(m,n); 
  iflip = false ;
  rflip = false ;
  for icol = 1:n
    mleft(:,icol) = icol + [0:n:n*n]' + (o-1)*sheet ;
  end
  for irow = start:(m-1-dr)
    for iblock=1:(n-1)
      esquare = zeros(2,2);
      esquare(1,1) = mleft(irow,iblock);
      esquare(1,2) = mleft(irow,iblock + 1);
      esquare(2,1) = mleft(irow + 1, iblock);
      esquare(2,2) = mleft(irow + 1, iblock + 1);
      if iflip
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,2);
        econn(enum,3) = esquare(1,1);
        econn(enum,4) = esquare(2,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,1);
        econn(enum,3) = esquare(2,2);
        econn(enum,4) = esquare(1,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      else
        econn(enum,1) = enum;
        %element 1
        econn(enum,2) = esquare(2,1);
        econn(enum,3) = esquare(1,2);
        econn(enum,4) = esquare(1,1);
        enum = enum + 1;
        %element 2
        econn(enum,1) = enum;
        econn(enum,2) = esquare(1,2);
        econn(enum,3) = esquare(2,1);
        econn(enum,4) = esquare(2,2);
        enum = enum + 1;
        %flip the block
        iflip = not(iflip);
      end
    end
    rflip = not(rflip);
    iflip = rflip;
  end
  econn(:,1) = econn(:,1) + 5*(m-1)*(n-1)*(o-1)*ones(nElem611,1);
  fprintf(fid,' %d   %d   %d   %d\n',econn'); 
end
fprintf(fid,'\n');
fprintf(fid,'\n');
if elType==305
% Boundary conditions
fprintf(fid,'%% BOUNDARY CONDITIONS\n');
fprintf(fid,'%%  In this section the number of lines is equal to one plus twice the number of nodes on\n');
fprintf(fid,'%%  the boundaries of the computational domain times the number of measurement sets. All\n');
fprintf(fid,'%%  the nodal dirichlet boundary conditions for one measurement set are given together\n');
fprintf(fid,'%%  before another measurement set is written.\n');
fprintf(fid,'%%  For every measurement set, the first line is the number of nodes on the boundary.\n');
fprintf(fid,'%%  The following lines work by pairs: the first line provides the node number and the\n');
fprintf(fid,'%%  second line indicates the type of boundary condition.\n');
fprintf(fid,'%%  On every second line, the numbers work by pairs: first is an integer set to 0 if the\n');
fprintf(fid,'%%  Dirichlet boundary is inactive and set to 1 if it is active. Inactive means that a\n');
fprintf(fid,'%%  homogeneous Neuman condition will be use instead of enforcing a displacement with a\n');
fprintf(fid,'%%  Dirichlet boundary. The second number of the pair is the value of the enforced\n');
fprintf(fid,'%%  displacement through the Dirichlet boundary condition.\n');
fprintf(fid,'%%  Here, the first pair is the displacement along the lateral direction, the second\n');
fprintf(fid,'%%  pair is the displacement along the axial direction and, when there is a third pair,\n');
fprintf(fid,'%%  it corresponds to a measure of the pressure inside the material.\n');
fprintf(fid,'boundaries\n');
boun=2*o*n+2*o*(m-2)+2*(n-2)*(m-2); % number of points belonging to the boundary
% define index numbers
tmp1=ones(m-2,1)*[2:1:(n-1)];% for face in elevational direction (z-axis)
tmp2=([2:1:(m-1)]')*ones(1,n-2);% for face in elevational direction (z-axis)
tmp3=ones(o,1)*[2:1:(m-1)];% for face in lateral direction (x-axis)
tmp4=([1:1:o]')*ones(1,m-2);% for face in lateral direction (x-axis)
tmp5=ones(o,1)*[1:1:n];% for face in axial direction (y-axis)
tmp6=([1:1:o]')*ones(1,n);% for face in axial direction (y-axis)
% top (close to probe), bottom (far from probe) and then sides (lateral and then elevational)
indexBoundi=[tmp5(:);tmp5(:);ones(o*(m-2),1);n*ones(o*(m-2),1);tmp1(:);tmp1(:)];
indexBoundj=[ones(o*n,1);m*ones(o*n,1);tmp3(:);tmp3(:);tmp2(:);tmp2(:)];
indexBoundk=[tmp6(:);tmp6(:);tmp4(:);tmp4(:);ones((n-2)*(m-2),1);o*ones((n-2)*(m-2),1)];
%indexBoundi=[ones(o*n,1);m*ones(o*n,1);tmp3(:);tmp3(:);tmp1(:);tmp1(:)];
%indexBoundj=[tmp5(:);tmp5(:);ones(o*(m-2),1);n*ones(o*(m-2),1);tmp2(:);tmp2(:)];
%indexBoundk=[tmp6(:);tmp6(:);tmp4(:);tmp4(:);ones((n-2)*(m-2),1);o*ones((n-2)*(m-2),1)];
inde=[indexBoundi+(indexBoundj-1)*n+(indexBoundk-1)*n*m];
% transpose nodesNum (to ndgrid format...)
nodesNumB=nodesNum(inde);
%inde=[indexBoundj+(indexBoundi-1)*m+(indexBoundk-1)*n*m];
%xCoor_bou=xCoor_mod(inde);% x coordinnate of the nodes on the boundary
%yCoor_bou=yCoor_mod(inde);% y coordinnate
indexAxialFaces=[1:1:(2*o*n)]';
indexLateralFaces=[(1+2*o*n):1:(2*o*n+2*o*(m-2))];
indexElevationalFaces=[(2*o*n+2*o*(m-2)+1):1:boun]';
indexElevationalLayers=cell(o,1);
%indexElevationalLayers{1,1}=[(2*o*n+2*o*(m-2)+1):1:(2*o*n+2*o*(m-2)+(n-2)*(m-2))];
for i=1:1:o
  indexElevationalLayers{i,1}=find(indexBoundk==i);
end
%indexElevationalLayers{o,1}=[(2*o*n+2*o*(m-2)+(n-2)*(m-2)+1):1:boun];
for lln=1:1:size(Ln,2)
  fprintf(fid,' %d\n',boun);
  %inde=[indexBoundj+(indexBoundi-1)*m+(indexBoundk-1)*m*n+(lln-1)*m*n*o];
  if elType==305
    % by default impose the axial displacement with DBC and NBC everywhere else
    inde=[indexBoundj+(indexBoundi-1)*m+(indexBoundk-1)*n*m+(lln-1)*o*m*n];
    bound_x_y_z=[nodesNumB,0*ones(boun,1),xDisp_mod(inde),1*ones(boun,1),yDisp_mod(inde),...
                 0*ones(boun,1),zDisp_mod(inde),0*ones(boun,1),pDisp_mod(inde)];
    if geoCylDispCart
      if false
        % smooth the displacements at the boundaries
        disp('BCs: x, y and z displacements are smoothed')
        xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
        xDisp_avg=create_input_smoothDisp(xDisp_avg);% smooth the displacements on the surface
        xDisp_avg=create_input_smoothDisp(xDisp_avg);% smooth the displacements on the surface
        bound_x_y_z(:,3)=xDisp_avg(inde);
        yDisp_avg=create_input_smoothDisp(yDisp_mod(:,:,:,lln));% smooth the displacements on the surface
        yDisp_avg=create_input_smoothDisp(yDisp_avg);% smooth the displacements on the surface
        yDisp_avg=create_input_smoothDisp(yDisp_avg);% smooth the displacements on the surface
        bound_x_y_z(:,5)=yDisp_avg(inde);
        zDisp_avg=create_input_smoothDisp(zDisp_mod(:,:,:,lln));% smooth the displacements on the surface
        zDisp_avg=create_input_smoothDisp(zDisp_avg);% smooth the displacements on the surface
        zDisp_avg=create_input_smoothDisp(zDisp_avg);% smooth the displacements on the surface
        bound_x_y_z(:,7)=zDisp_avg(inde);
      else
        % recover the original axial displacement and use this one on the y-axis
	disp('BCs: recovering axial, lateral and elevational instead of cartesian');
        for k=1:1:size(xCoor_mod,3)
          % build the weighting matrix
          lateral=[xCoor_mod(2,1,k)-xCoor_mod(1,1,k),...
                   yCoor_mod(2,1,k)-yCoor_mod(1,1,k),...
                   zCoor_mod(2,1,k)-zCoor_mod(1,1,k)];
          if lateral(1)<0
            lateral=-lateral;
          end
          lateral=lateral/sqrt(sum(lateral.^2,2));
          axial=[xCoor_mod(1,2,k)-xCoor_mod(1,1,k),...
                 yCoor_mod(1,2,k)-yCoor_mod(1,1,k),...
                 zCoor_mod(1,2,k)-zCoor_mod(1,1,k)];
          if axial(2)<0
            axial=-axial;
          end
          axial=axial/sqrt(sum(axial.^2,2));
          elevational=cross(lateral,axial);
          elevational=elevational/sqrt(sum(elevational.^2,2));
          if enforceOrthogonality
	    lateral=cross(axial,elevational);
            lateral=lateral/sqrt(sum(lateral.^2,2));
	  end
	  axVect(k,1:3)=axial;
	  laVect(k,1:3)=lateral;
	  elVect(k,1:3)=elevational;
	end
	xDisp_avg=zeros(size(yDisp_mod,1),size(yDisp_mod,2),size(yDisp_mod,3));
        for k=1:1:size(yDisp_mod,3)
	  selectMat=(laVect(k,1:3)')*laVect(k,1:3);
	  for j=1:1:size(yDisp_mod,1)
	    for i=1:1:size(yDisp_mod,2)
	      projVect=[xDisp_mod(j,i,k,lln),yDisp_mod(j,i,k,lln),zDisp_mod(j,i,k,lln)];
	      vect=selectMat*(projVect');
              xDisp_avg(j,i,k)=sqrt(sum(vect.^2,1));
	      if vect(1)<0
	        xDisp_avg(j,i,k)=-xDisp_avg(j,i,k);
	      end
	    end
	  end
	end
	yDisp_avg=zeros(size(yDisp_mod,1),size(yDisp_mod,2),size(yDisp_mod,3));
        for k=1:1:size(yDisp_mod,3)
	  selectMat=(axVect(k,1:3)')*axVect(k,1:3);
	  for j=1:1:size(yDisp_mod,1)
	    for i=1:1:size(yDisp_mod,2)
	      projVect=[xDisp_mod(j,i,k,lln),yDisp_mod(j,i,k,lln),zDisp_mod(j,i,k,lln)];
	      vect=selectMat*(projVect');
              yDisp_avg(j,i,k)=sqrt(sum(vect.^2,1));
	      if vect(2)<0
	        yDisp_avg(j,i,k)=-yDisp_avg(j,i,k);
	      end
	    end
	  end
	end
	zDisp_avg=zeros(size(yDisp_mod,1),size(yDisp_mod,2),size(yDisp_mod,3));
        for k=1:1:size(yDisp_mod,3)
	  selectMat=(elVect(k,1:3)')*elVect(k,1:3);
	  for j=1:1:size(yDisp_mod,1)
	    for i=1:1:size(yDisp_mod,2)
	      projVect=[xDisp_mod(j,i,k,lln),yDisp_mod(j,i,k,lln),zDisp_mod(j,i,k,lln)];
	      vect=selectMat*(projVect');
              zDisp_avg(j,i,k)=sqrt(sum(vect.^2,1));
	      if vect(3)<0
	        zDisp_avg(j,i,k)=-zDisp_avg(j,i,k);
	      end
	    end
	  end
	end
        bound_x_y_z(:,3)=xDisp_avg(inde);
        bound_x_y_z(:,5)=yDisp_avg(inde);
        bound_x_y_z(:,7)=zDisp_avg(inde);
	clear laVect; clear axVect; clear elVect; clear selectMat; clear projVect; clear vect;
      end
    end
    inde=[indexBoundj+(indexBoundi-1)*m+(indexBoundk-1)*n*m];
    %input('here')
    if true% impose DBC in axial, homogeneous NBC for lateral & elevational (with 3 Dofs=DBC against RBM)
      disp('BCs: imposing NBC for lateral and elevational');
      % to avoid rigid body modes, look for 2 points close to the center of the face indexBoundj==1
      %  where the lateral displacement is close to the average of the lateral disp; repeat for elevational
      %if (m>5.5) mm=max([1,floor(m/6)],[],2); else mm=0; end% not needed 
      if (n>5.5) nn=max([1,floor(n/6)],[],2); else nn=0; end
      if (o>5.5) oo=max([1,floor(o/6)],[],2); else oo=0; end
      indexBoundi2=tmp5([(1+oo):(o-oo)],[(1+nn):(n-nn)]);
      indexBoundj2=ones(o-2*oo,n-2*nn);
      indexBoundk2=tmp6([(1+oo):(o-oo)],[(1+nn):(n-nn)]);
      inds3=indexBoundj2+(indexBoundi2-1)*m+(indexBoundk2-1)*m*n;
      inds=inds3(:);
      xDisp_sel=zeros(size(inds3));% selected window
      xDisp_sel(:)=xDisp_mod(inds+(lln-1)*n*m*o);
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      xDisp_avg_sel=zeros(size(inds3));
      xDisp_avg_sel(:)=xDisp_avg(inds);
      tmp1=abs(xDisp_avg_sel-xDisp_sel);
      zDisp_sel=zeros(size(inds3));% selected window
      zDisp_sel(:)=zDisp_mod(inds+(lln-1)*n*m*o);
      zDisp_avg=create_input_smoothDisp(zDisp_mod(:,:,:,lln));% average field
      zDisp_avg_sel=zeros(size(inds3));
      zDisp_avg_sel(:)=zDisp_avg(inds);
      tmp2=abs(zDisp_avg_sel-zDisp_sel);
      disp('BCs:  fixing 3 DOFs against RBM');
      tmp=tmp1(:)+tmp2(:);
      [val,inds2]=sort(tmp,'ascend');
      inds2=inds2(1);% les deux points sont-ils assez espaces?
      ind5=find(inde==inds(inds2));
      bound_x_y_z(ind5,2)=1;
      %bound_x_y_z(ind5,3)=xDisp_avg_sel(inds2);
      bound_x_y_z(ind5,6)=1;
      %bound_x_y_z(ind5,7)=zDisp_avg_sel(inds2);
      % find the best point in x direction and fix the zDisp
      tmp=floor((inds2-1)/(o-2*oo));
      inds6=[(inds2-tmp*(o-2*oo)):(o-2*oo):(inds2-tmp*(o-2*oo)+(o-2*oo)*(n-2*nn-1))];
      inds6=setdiff(inds6,(inds2-3*(o-2*oo)));% take a point far enough
      inds6=setdiff(inds6,(inds2-2*(o-2*oo)));
      inds6=setdiff(inds6,(inds2-(o-2*oo)));
      inds6=setdiff(inds6,inds2);
      inds6=setdiff(inds6,(inds2+(o-2*oo)));
      inds6=setdiff(inds6,(inds2+2*(o-2*oo)));
      inds6=setdiff(inds6,(inds2+3*(o-2*oo)));% if inds6 is empty here, the domain is too small
      tmp3=tmp2(inds6);
      tmp=tmp3(:);
      [val,inds2]=sort(tmp,'ascend');
      ind5=find(inde==inds(inds6(inds2(1))));
      bound_x_y_z(ind5,6)=1;
      %bound_x_y_z(ind5,7)=zDisp_avg_sel(inds4);
    elseif false
      % fix DBC=-0.268 everywhere for elevational disp and 1 DBC for axial to avoid RBM
      %if (m>5.5) mm=max([1,floor(m/6)],[],2); else mm=0; end% not needed 
      if (n>5.5) nn=max([1,floor(n/6)],[],2); else nn=0; end
      if (o>5.5) oo=max([1,floor(o/6)],[],2); else oo=0; end
      inds3=ones(o-2*oo,n-2*nn)+(tmp5([(1+oo):(o-oo)],[(1+nn):(n-nn)])-1)*m+(tmp6([(1+oo):(o-oo)],[(1+nn):(n-nn)])-1)*m*n;
      inds=inds3(:);
      xDisp_sel=zeros(size(inds3));% selected window
      xDisp_sel(:)=xDisp_mod(inds+(lln-1)*n*m*o);
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      xDisp_avg_sel=zeros(size(inds3));
      xDisp_avg_sel(:)=xDisp_avg(inds);
      tmp=abs(xDisp_avg_sel-xDisp_sel);
      tmp=tmp(:);
      [val,inds2]=sort(tmp,'ascend');
      %inds4=inds2(2);
      inds2=inds2(1);% les deux points sont-ils assez espaces?
      ind5=find(inde==inds(inds2));
      bound_x_y_z(ind5,2)=1;
      %ind5=find(inde==inds(inds4));
      %bound_x_y_z(ind5,2)=1;
      bound_x_y_z(:,6)=1;
      bound_x_y_z(:,7)=-0.268;
    elseif false
      % fix DBC=-0.268 on the elevational faces for elevational disp and 1 DBC for axial to avoid RBM
      %if (m>5.5) mm=max([1,floor(m/6)],[],2); else mm=0; end% not needed 
      if (n>5.5) nn=max([1,floor(n/6)],[],2); else nn=0; end
      if (o>5.5) oo=max([1,floor(o/6)],[],2); else oo=0; end
      inds3=ones(o-2*oo,n-2*nn)+(tmp5([(1+oo):(o-oo)],[(1+nn):(n-nn)])-1)*m+(tmp6([(1+oo):(o-oo)],[(1+nn):(n-nn)])-1)*m*n;
      inds=inds3(:);
      xDisp_sel=zeros(size(inds3));% selected window
      xDisp_sel(:)=xDisp_mod(inds+(lln-1)*n*m*o);
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      xDisp_avg_sel=zeros(size(inds3));
      xDisp_avg_sel(:)=xDisp_avg(inds);
      tmp=abs(xDisp_avg_sel-xDisp_sel);
      tmp=tmp(:);
      [val,inds2]=sort(tmp,'ascend');
      %inds4=inds2(2);
      inds2=inds2(1);% les deux points sont-ils assez espaces?
      ind5=find(inde==inds(inds2));
      bound_x_y_z(ind5,2)=1;
      %ind5=find(inde==inds(inds4));
      %bound_x_y_z(ind5,2)=1;
      bound_x_y_z(indexElevationalFaces,6)=1;
      bound_x_y_z(:,7)=-0.268;
    elseif false% assad1
      disp('BCs: DBC=(mean of elev. disp on elev faces) for elevational dofs everywhere');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on axial faces');
      bound_x_y_z(:,6)=1;
      tmp=zDisp_mod(:,:,1,lln);
      tmp2=zDisp_mod(:,:,size(zDisp_mod,3),lln);
      bound_x_y_z(:,7)=(mean(tmp(:))+mean(tmp2(:)))/2;
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
     %elseif true% assad12
    elseif false %tom
      disp('BCs: DBC=smoothed measured disp for elevational dofs everywhere');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on axial faces');
      bound_x_y_z(:,6)=1;
      zDisp_avg=create_input_smoothDisp(zDisp_mod(:,:,:,lln));
      bound_x_y_z(:,7)=zDisp_avg(inde);
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false% assad13
      disp('BCs: DBC=smoothed measured disp for elevational dofs on axial+elev faces');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on axial+lateral faces');
      zDisp_avg=create_input_smoothDisp(zDisp_mod(:,:,:,lln));
      bound_x_y_z(:,7)=zDisp_avg(inde);
      bound_x_y_z(indexAxialFaces,6)=1;
      bound_x_y_z(indexElevationalFaces,6)=1;
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(:,3)=xDisp_avg(inde);
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(indexLateralFaces,2)=1;
    elseif false% assad2
      disp('BCs: DBC=(mean of elev. disp on elev faces) for elevational dofs everywhere');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on elevational faces');
      bound_x_y_z(:,6)=1;
      tmp=zDisp_mod(:,:,1,lln);
      tmp2=zDisp_mod(:,:,size(zDisp_mod,3),lln);
      bound_x_y_z(:,7)=(mean(tmp(:))+mean(tmp2(:)))/2;
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexElevationalFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif true% paul1
      disp('BCs: DBC=measured disp for elevational dofs everywhere');
      disp('BCs: DBC=measured disp for lateral dofs on axial+elevational faces');
      bound_x_y_z(:,6)=1;
      zDisp_avg=zDisp_mod(:,:,:,lln);%create_input_smoothDisp(zDisp_mod(:,:,:,lln));
      bound_x_y_z(:,7)=zDisp_avg(inde);
      xDisp_avg=xDisp_mod(:,:,:,lln);%create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(indexElevationalFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false% paul2
      disp('BCs: DBC=measured disp for elevational dofs on axial+elevational faces');
      disp('BCs: DBC=measured disp for lateral dofs on axial faces');
      bound_x_y_z(indexAxialFaces,6)=1;
      bound_x_y_z(indexElevationalFaces,6)=1;
      zDisp_avg=zDisp_mod(:,:,:,lln);%create_input_smoothDisp(zDisp_mod(:,:,:,lln));
      bound_x_y_z(:,7)=zDisp_avg(inde);
      xDisp_avg=xDisp_mod(:,:,:,lln);%create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false% paul3
      disp('BCs: DBC=measured disp for elevational dofs on axial+elevational faces');
      disp('BCs: DBC=measured disp for lateral dofs on axial faces');
      bound_x_y_z(indexAxialFaces,6)=1;
      bound_x_y_z(indexElevationalFaces,6)=1;
      zDisp_avg=zDisp_mod(:,:,:,lln);%create_input_smoothDisp(zDisp_mod(:,:,:,lln));
      bound_x_y_z(:,7)=zDisp_avg(inde);
      xDisp_avg=xDisp_mod(:,:,:,lln);%create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(indexElevationalFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false% JF0
      disp('BCs: DBC=(mean of elev. disp on each elev face) for elevational dofs on each elev face and linear interpolation otherwise');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on axial faces');
      bound_x_y_z(:,6)=1;
      tmp=zDisp_mod(:,:,1,lln);
      tmp2=zDisp_mod(:,:,size(zDisp_mod,3),lln);
      for i=1:1:o
        bound_x_y_z(indexElevationalLayers{i,1},7)=((o-i)/(o-1))*mean(tmp(:))+((i-1)/(o-1))*mean(tmp2(:));
      end
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false% JF1
      disp('BCs: DBC=(mean of elev. disp on elev faces) for elevational dofs everywhere');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on axial+elevational faces');
      bound_x_y_z(:,6)=1;
      tmp=zDisp_mod(:,:,1,lln);
      tmp2=zDisp_mod(:,:,size(zDisp_mod,3),lln);
      bound_x_y_z(:,7)=(mean(tmp(:))+mean(tmp2(:)))/2;
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(indexElevationalFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false% JF2
      disp('BCs: DBC=(mean of elev. disp on each elev face) for elevational dofs on each elev face and linear interpolation otherwise');
      disp('BCs: DBC=smoothed measured disp for lateral dofs on axial+elevational faces');
      bound_x_y_z(:,6)=1;
      tmp=zDisp_mod(:,:,1,lln);
      tmp2=zDisp_mod(:,:,size(zDisp_mod,3),lln);
      for i=1:1:o
        bound_x_y_z(indexElevationalLayers{i,1},7)=((o-i)/(o-1))*mean(tmp(:))+((i-1)/(o-1))*mean(tmp2(:));
      end
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(indexAxialFaces,2)=1;
      bound_x_y_z(indexElevationalFaces,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
    elseif false
      disp('BCs: DBC=(mean of elev. disp on elev faces) for elevational dofs everywhere');
      disp('BCs: DBC=smoothed measued disp for lateral dofs on all faces');
      bound_x_y_z(:,6)=1;
      tmp=zDisp_mod(:,:,1,lln);
      tmp2=zDisp_mod(:,:,size(zDisp_mod,3),lln);
      bound_x_y_z(:,7)=(mean(tmp(:))+mean(tmp2(:)))/2;
      xDisp_avg=create_input_smoothDisp(xDisp_mod(:,:,:,lln));% smooth the displacements on the surface
      bound_x_y_z(:,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
      % plus impose a pressure DOF with DBC if all other DOFs are constrained
      bound_x_y_z(1,8)=1;
    else
      % impose DBC everywhere using smoothed lateral and elevational displacements
      disp('BCs: impose DBC for every DOF with smooth measured displacements');
      xDisp_avg2=create_input_smoothDisp(xDisp_mod(:,:,:,lln));
      xDisp_avg3=create_input_smoothDisp(xDisp_avg2);
      xDisp_avg=create_input_smoothDisp(xDisp_avg3);
      bound_x_y_z(:,2)=1;
      bound_x_y_z(:,3)=xDisp_avg(inde);
      zDisp_avg2=create_input_smoothDisp(zDisp_mod(:,:,:,lln));
      zDisp_avg3=create_input_smoothDisp(zDisp_avg2);
      zDisp_avg4=create_input_smoothDisp(zDisp_avg3);
      zDisp_avg5=create_input_smoothDisp(zDisp_avg4);
      zDisp_avg=create_input_smoothDisp(zDisp_avg5);
      bound_x_y_z(:,6)=1;
      bound_x_y_z(:,7)=zDisp_avg(inde);
      % plus impose a pressure DOF with DBC if all other DOFs are constrained
      bound_x_y_z(1,8)=1;
    end
    fprintf(fid, '  %d\n   %d %11.6f   %d %11.6f   %d %11.6f  %d %11.6f\n',bound_x_y_z');
  end
end
end% if elType==305


if elType==305
  % Data initial value
  fprintf(fid,'\n');
  fprintf(fid,'\n');
  fprintf(fid,'%% NODAL VALUE OF A FIELD\n');
  fprintf(fid,'%%  This section has as many lines as there are nodes plus one; the order for the\n');
  fprintf(fid,'%%  nodal data is the same as in section ''coor''.\n');
  fprintf(fid,'%%  The first line is the number of parameters at every node.\n');
  fprintf(fid,'%%  On every subsequent line, the numbers work by pairs: first is an integer set to\n');
  fprintf(fid,'%%  0 if the parameter is inactive (not optimized and kept to its initial value) or\n');
  fprintf(fid,'%%  set to 1 if the parameter is an optimization variable. The second number of the\n');
  fprintf(fid,'%%  pair is the initial value of the parameter.\n');
  fprintf(fid,'%%  Here, these parameters describe the stress-strain relation in the material, the\n');
  fprintf(fid,'%%  first parameter is gamma (non-linearity) and the second is mu (shear modulus) of\n');
  fprintf(fid,'%%  the Veronda-Westman model.\n');
  fprintf(fid,'datn\n');
  uu=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
  vv=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
  ww=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
  xx=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
  if elType==305
    real_disp=zeros(max_node,8);% matrix to store the printed field
  elseif elType==31||elType==32
    real_disp=zeros(max_node,6);% matrix to store the printed field
  end
  for lln=1:1:size(Ln,2)
    for k=1:1:size(xDisp_mod,3)
      uu(:,:,k)=xDisp_mod(:,:,k,lln)';
      vv(:,:,k)=yDisp_mod(:,:,k,lln)';
      ww(:,:,k)=zDisp_mod(:,:,k,lln)';
      xx(:,:,k)=pDisp_mod(:,:,k,lln)';
    end
  end
  if not(addElem611)
    datn=2;
    fprintf(fid,' %d\n',datn);
    data(1:max_node,1)=0;% active/inactive optimization parameter 1
    data(1:max_node,2)=0.0;% initial value parameter 1
    data(1:max_node,3)=1;% active/inactive optimization parameter 2
    data(1:max_node,4)=1.05;% initial value parameter 2
      fprintf(fid,'  %d %11.6f   %d %11.6f\n',data');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
  else
    % %datn=2;
    % datn=4;% make 5 if nonlinear;
    % fprintf(fid,' %d\n',datn);
    % %data(1:max_node,1)=0;% active/inactive optimization parameter 1
    % %data(1:max_node,2)=0.0;% initial value parameter 1
    % data(1:max_node,1)=1;% active/inactive optimization parameter 2
    % data(1:max_node,2)=1.05;% initial value parameter 2
    % data(1:max_node,3)=0;
    % data(1:max_node,4)=uu(:);
    % data(1:max_node,5)=0;
    % data(1:max_node,6)=vv(:);
    % data(1:max_node,7)=0;
    % data(1:max_node,8)=ww(:);
  
    %nonlinear options
    datn=5;% make 5 if nonlinear;
    fprintf(fid,' %d\n',datn);
    data(1:max_node,1)=0;% active/inactive optimization parameter 1
    data(1:max_node,2)=0.0;% initial value parameter 1
    data(1:max_node,3)=1;% active/inactive optimization parameter 2
    data(1:max_node,4)=1.05;% initial value parameter 2
    data(1:max_node,5)=0;
    data(1:max_node,6)=uu(:);
    data(1:max_node,7)=0;
    data(1:max_node,8)=vv(:);
    data(1:max_node,9)=0;
    data(1:max_node,10)=ww(:);
    fprintf(fid,'  %d %11.6f   %d %11.6f   %d %11.6f   %d %11.6f   %d %11.6f\n',data');
    %fprintf(fid,'  %d %11.6f   %d %11.6f   %d %11.6f   %d %11.6f\n',data');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
  end
elseif elType==32
  % enter the rotation matrix for each point (uniform over each elevational plane when appropriate)
  rotations=zeros(3,3,size(xCoor_mod,3));
  for k=1:1:size(xCoor_mod,3)
    lateral=[xCoor_mod(2,1,k)-xCoor_mod(1,1,k);...
           yCoor_mod(2,1,k)-yCoor_mod(1,1,k);...
	   zCoor_mod(2,1,k)-zCoor_mod(1,1,k)];
    lateral=lateral/sqrt(sum(lateral.^2,1));
    if lateral(1)<0
      lateral=-lateral;
    end
    axial=[xCoor_mod(1,2,k)-xCoor_mod(1,1,k);...
           yCoor_mod(1,2,k)-yCoor_mod(1,1,k);...
	   zCoor_mod(1,2,k)-zCoor_mod(1,1,k)];
    axial=axial/sqrt(sum(axial.^2,1));
    if axial(2)<0
      axial=-axial;
    end
    elev=cross(lateral,axial);
    elev=elev/sqrt(sum(elev.^2,1));
%k
    MM=[lateral,axial,elev];
    %MM=(MM-MM')/2;
%MM
%[MM*(MM'),(MM')*MM]
%MM^-1
    rotations(:,:,k)=MM;% rotation from cartesian to cylindrical coordinates
  end
  MM11=zeros(size(xCoor_mod));
  MM12=zeros(size(xCoor_mod));
  MM13=zeros(size(xCoor_mod));
  MM21=zeros(size(xCoor_mod));
  MM22=zeros(size(xCoor_mod));
  MM23=zeros(size(xCoor_mod));
  MM31=zeros(size(xCoor_mod));
  MM32=zeros(size(xCoor_mod));
  MM33=zeros(size(xCoor_mod));
  for k=1:1:size(xCoor_mod,3)
    MM11(:,:,k)=rotations(1,1,k);
    MM12(:,:,k)=rotations(1,2,k);
    MM13(:,:,k)=rotations(1,3,k);
    MM21(:,:,k)=rotations(2,1,k);
    MM22(:,:,k)=rotations(2,2,k);
    MM23(:,:,k)=rotations(2,3,k);
    MM31(:,:,k)=rotations(3,1,k);
    MM32(:,:,k)=rotations(3,2,k);
    MM33(:,:,k)=rotations(3,3,k);
  end
  on=ones(size(MM33(:)));
  MMM=[on,MM11(:),on,MM12(:),on,MM13(:),on,MM21(:),on,MM22(:),on,MM23(:),on,MM31(:),on,MM32(:),on,MM33(:)];
  fprintf(fid,'%% NODAL ROTATION MATRIX\n');
  fprintf(fid,'%%  to match data given in non cartesian coordinates\n');
  fprintf(fid,'datn\n');
  datn=9;
  fprintf(fid,' %d\n',datn);
  fprintf(fid,'%d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f %d %f\n',MMM');
fprintf(fid,'\n');
fprintf(fid,'\n');
else
  disp('such an elType is not implemented for te datn section ... exiting');
  return
end% if elType==305


fprintf(fid,'%% MEASURED DISPLACEMENTS\n');
fprintf(fid,'%%  This section has a number of lines equal to the number of nodes times the number of\n');
fprintf(fid,'%%  measurement sets used for the optimization. All the nodal displacements for one\n');
fprintf(fid,'%%  measurement set are given together before another measurement set is written. Within\n');
fprintf(fid,'%%  every measurement set, the displacements appear in the same order as the nodes in\n');
fprintf(fid,'%%  section ''coor''.\n');
fprintf(fid,'%%  On every line, the numbers work by pairs: first is a real number that is a weight\n');
fprintf(fid,'%%  for the associated displacement when computing the objective function. The second\n');
fprintf(fid,'%%  number of the pair is the measured displacement.\n');
fprintf(fid,'%%  Here, the first pair is the displacement along the lateral direction, the second\n');
fprintf(fid,'%%  pair is the displacement along the axial direction and, when there is a fourth pair,\n');
fprintf(fid,'%%  it corresponds to a measure of the pressure inside the material.\n');
fprintf(fid,'measuredDisplacements\n');
uu=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
vv=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
ww=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
xx=zeros(size(xDisp_mod,2),size(xDisp_mod,1),size(xDisp_mod,3));% intermediate matrix
if elType==305
  real_disp=zeros(max_node,8);% matrix to store the printed field
elseif elType==31||elType==32
  real_disp=zeros(max_node,6);% matrix to store the printed field
end
for lln=1:1:size(Ln,2)
  for k=1:1:size(xDisp_mod,3)
    uu(:,:,k)=xDisp_mod(:,:,k,lln)';
    vv(:,:,k)=yDisp_mod(:,:,k,lln)';
    ww(:,:,k)=zDisp_mod(:,:,k,lln)';
    xx(:,:,k)=pDisp_mod(:,:,k,lln)';
  end

  if elType==305
    % lateral displacement
    if false% weight of the lateral displacement in the objective function
      % make a variable map depending on teh position of the point
      disp('the weight of the lateral disp in the obj func depends on the position of the node');
      si=size(xDisp_mod);
      nRemove=6;% layers for every boundary (top, bottom,side)
      nTransition=4;
      maxWeight=0.3;
      tmp=zeros(si(1),si(2));
      for ii=1:1:nTransition
        tmp((1+nRemove+ii):(si(1)-nRemove-ii),(1+nRemove+ii):(si(2)-nRemove-ii))=maxWeight*ii/nTransition;
      end
      real_disp(:,1)=tmp(:);
    else
      % set everywhere the same value
      real_disp(:,1)=0.0;
    end
    % axial displacement
    if lln==1
      real_disp(:,3)=1.0;
    else
      real_disp(:,3)=1.0;
    end
    % elevational displacement
    real_disp(:,5)=0.0;
    % pressure displacement
    real_disp(:,7)=0.0;
    if geoCylDispCart
      wAxial=1.0;
      wLateral=0.0;
      wElevational=0.0;
      for k=1:1:size(xCoor_mod,3)
        % build the weighting matrix
        lateral=[xCoor_mod(2,1,k)-xCoor_mod(1,1,k),...
                 yCoor_mod(2,1,k)-yCoor_mod(1,1,k),...
                 zCoor_mod(2,1,k)-zCoor_mod(1,1,k)];
        if lateral(1)<0
          lateral=-lateral;
        end
        lateral=lateral/sqrt(sum(lateral.^2,2));
        axial=[xCoor_mod(1,2,k)-xCoor_mod(1,1,k),...
               yCoor_mod(1,2,k)-yCoor_mod(1,1,k),...
               zCoor_mod(1,2,k)-zCoor_mod(1,1,k)];
        if axial(2)<0
          axial=-axial;
        end
        axial=axial/sqrt(sum(axial.^2,2));
        elevational=cross(lateral,axial);
        elevational=elevational/sqrt(sum(elevational.^2,2));
	if enforceOrthogonality
	  lateral=cross(axial,elevational);
          lateral=lateral/sqrt(sum(lateral.^2,2));
	end
	wMat=wLateral*((lateral')*lateral)+...
	     wAxial*((axial')*axial)+...
	     wElevational*((elevational')*elevational);
        % select the diagonal terms
	wDiag(1,k)=wMat(1,1);
	wDiag(2,k)=wMat(2,2);
	wDiag(3,k)=wMat(3,3);
	real_disp([1:m*n]+(k-1)*m*n,1)=wMat(1,1);
	real_disp([1:m*n]+(k-1)*m*n,3)=wMat(2,2);
	real_disp([1:m*n]+(k-1)*m*n,5)=wMat(3,3);
        % store the off-diagonal terms
	wOffD(1,k)=wMat(1,2);
	wOffD(2,k)=wMat(1,3);
	wOffD(3,k)=wMat(2,3);
      end
    end
  elseif elType==31||elType==32
    real_disp(:,1)=1.0;% should not matter
    real_disp(:,3)=1.0;% should not matter
    real_disp(:,5)=1.0;% should not matter
  end


  real_disp(:,2)=uu(:);
  real_disp(:,4)=vv(:);
  real_disp(:,6)=ww(:);
  if elType==305
    real_disp(:,8)=xx(:);
  end
  if elType==305
    fprintf(fid,'%10.5f %11.6f %10.5f %11.6f %10.5f %11.6f %10.5f %11.6f\n', real_disp');
  elseif elType==31||elType==32
    fprintf(fid,'%8.3f %11.6f %8.3f %11.6f %8.3f %11.6f\n', real_disp');
  else%if elType==606
    disp('pb here ... exiting');
    return
    %fprintf(fid,'%8.3f %11.6f %8.3f %11.6f\n', real_disp(:,[1:4])');
  end
end
fprintf(fid,'\n');
fprintf(fid,'\n');

if geoCylDispCart% print the off-diagonal entries of the weight matrix
  for k=1:1:size(xCoor_mod,3)
    real_disp([1:1:m*n]+(k-1)*m*n,1)=wOffD(1,k);
    real_disp([1:1:m*n]+(k-1)*m*n,2)=wOffD(2,k);
    real_disp([1:1:m*n]+(k-1)*m*n,3)=0.0;% pressure dof
    real_disp([1:1:m*n]+(k-1)*m*n,4)=wOffD(3,k);
    real_disp([1:1:m*n]+(k-1)*m*n,5)=0.0;% pressure dof
    real_disp([1:1:m*n]+(k-1)*m*n,6)=0.0;% pressure dof
  end
  fprintf(fid,'%% Off-diagonal terms for the weighting of displacement when the\n');
  fprintf(fid,'%%  grid is structured but irregular\n');
  fprintf(fid,'fmnd\n');
  fprintf(fid,' %f %f %f %f %f %f\n',real_disp(:,[1:6])');
  fprintf(fid,'\n');
  fprintf(fid,'\n');
end

fprintf(fid,'end\n');

fclose(fid);

%if plotBCs
%xDispi=xDispi(ay,ax,Ln);
%xDispNlm=xDispNlm(ay,ax,Ln);
%xDispNlms=xDispNlms(ay,ax,Ln);
%for lln=1:1:size(Ln,2)
%  indexes=indexBoundi+(indexBoundj-1)*m+n*m*(lln-1);
%  boundi{lln}=xDispi(indexes);
%  boundNlm{lln}=xDispNlm(indexes);
%  boundNlms{lln}=xDispNlms(indexes);
%end
%%inds=[[1:1:n],[(2*n+m-2+1):1:(2*n+2*(m-2))],[(2*n):-1:(n+1)],[(2*n+m-2):-1:(2*n+1)]];
%for lln=1:1:size(Ln,2)
%figure
%set(gcf,'color','w');
%hold on
%plot(boundi{lln}(inds),'k-','linewidth',2);
%plot(boundNlm{lln}(inds),'b-','linewidth',2);
%plot(boundNlms{lln}(inds),'r-','linewidth',2);
%plot(boundx{lln}(inds),'g-','linewidth',2);
%legend('initial','nlmeans','nlmeans+laplace','linearized');
%title(['lateral displacements at the boundaries for frame ' num2str(Ln(1,lln))]);
%end
%end

