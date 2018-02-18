%clc;
close all;
clear all;

%Loading the grid coordinates and the measured displacement at the right format (cf README)
%filename='CC024summedDisplacements21_35f.mat';
filepath='';% default
filename='target1_compression_accumfsm.mat';
%filename='syntheticSevan606-noise6.mat';%'target4_compression_accum3f.mat';
%filename='3D22DdataThin.mat'
%filename='3D22DdataAll.mat'
%filename='3D22DdataThinGap.mat'
%filename='3D22DdataThinGapWithProbe.mat'
%filename='3D22Ddata.mat'
%filename='3D22Ddata45-2per.mat'
%filename='3D22Ddata45-5per.mat'
%filename='3D22DdataInvertedGamma1.mat'
%filename='3D22DdataNoContrastGamma.mat'
%filename='dseries_accum_110207_1f.mat'; filepath='../database3D/';
%filename='nonlinph-accum-20110218-ud3f.mat'; filepath='../database3D/';
%filename='forward606-step-noise5precent.mat';

% select a series of loading steps (also called frames)
Ln=[1];%[1];%[3,15,25,35];

% select every "reduce" measured point
reduce=5; %30 for tiny size; %1 for highest resolution
% if the data on its boundary is inaccurate, define margins (in number of lines/columns)
topMargin=40;% closest to the ultrasound transducer
bottomMargin=10;% far away from transducer
sideMargin=9;% left and right sides are treated equally

% choose a type of surfacic element
elType=66;%505, 606, 66 or 607

%add spring boundary elements
addElem611=true;
ksides=4; % put springs on 2,3, or 4 sides (right/left, bottom, then top)

elOrder=2;% for elem606 and elem608 only, elOrder=2 is linear element

plotBCs=false;% for debug or information

% if the data is 3D:
elevCut=true;% =true: make a slice with normal=elevationAxis(0,0,1); =false: normal=lateralAxis(0,1,0)
sliceNum=25;% plane number chosen to make slice


%----------------------- creating the input file ---------------%

% if the displacements have been filtered, you can load the unfiltered displacements here and save them
%  otherwise, just load filename...
if plotBCs
  disp('create_input: Check the filename hereafter...');
  %load([filepath 'syntheticSevan606-noise6.mat']);
  %load([filepath 'target1_compression_accumf.mat']);
  load([filepath 'target1_compression_accumfsm.mat']);
  %load 3D22Ddata.mat;
  %load 3D22Ddata454545.mat;
  xDispi=xDisp;  yDispi=yDisp;
end
disp(['create_input: Loading dataset: ' filepath filename]);
load([filepath filename]);
%yDisp(:,:,1)=yDisp(:,:,5)-yDisp(:,:,3);% JFD made for testing the quality of the data
%yDisp(:,:,2)=yDisp(:,:,7)-yDisp(:,:,3);
%xDisp(:,:,1)=xDisp(:,:,5)-xDisp(:,:,3);
%xDisp(:,:,2)=xDisp(:,:,7)-xDisp(:,:,3);

if nDim==3
  % reduce the data to match the 2D format
  disp(['create_input: This is 3D data; making a slice.']);
  if elevCut% sliceNormal=elevationAxis(0,0,1)
    xCoort=xCoor(:,:,sliceNum);
    yCoort=yCoor(:,:,sliceNum);
    xDispt=zeros(size(xDisp,1),size(xDisp,2),size(xDisp,4));
    yDispt=zeros(size(xDisp,1),size(xDisp,2),size(xDisp,4));
    pDispt=zeros(size(xDisp,1),size(xDisp,2),size(xDisp,4));
    for i=1:1:size(xDisp,4)
      xDispt(:,:,i)=xDisp(:,:,sliceNum,i);
      yDispt(:,:,i)=yDisp(:,:,sliceNum,i);
      pDispt(:,:,i)=pDisp(:,:,sliceNum,i);
    end
    xDisp=xDispt;
    yDisp=yDispt;
    xCoor=xCoort;
    yCoor=yCoort;
    pDisp=pDispt;
    clear xDispt; clear yDispt; clear pDispt; clear xCoort ;clear yCoort;
    clear zDisp; clear zCoor;
  else
    xCoort=zeros(size(xCoor,1),size(xCoor,3));
    aa=zCoor(:,sliceNum,:);
    xCoort(:)=aa(:);
    yCoort=zeros(size(yCoor,1),size(yCoor,3));
    aa=yCoor(:,sliceNum,:);
    yCoort(:)=aa(:);
    xDispt=zeros(size(xDisp,1),size(xDisp,3),size(xDisp,4));
    yDispt=zeros(size(xDisp,1),size(xDisp,3),size(xDisp,4));
    pDispt=zeros(size(xDisp,1),size(xDisp,3),size(xDisp,4));
    for i=1:1:size(xDisp,4)
      for j=1:1:size(xDisp,3)
        xDispt(:,j,i)=zDisp(:,sliceNum,j,i);
        yDispt(:,j,i)=yDisp(:,sliceNum,j,i);
        pDispt(:,j,i)=pDisp(:,sliceNum,j,i);
      end
    end
    xDisp=xDispt;
    yDisp=yDispt;
    xCoor=xCoort;
    yCoor=yCoort;
    pDisp=pDispt;
    clear xDispt; clear yDispt; clear pDispt; clear xCoort ;clear yCoort; clear aa;
    clear zDisp; clear zCoor;
  end
end

if plotBCs  xDispFiltered=xDisp; end

if addElem611
  disp(['create_input: Adding spring elements on ' num2str(ksides) ' sides']);
end

% size of the data
ex = size(xDisp,2);
ey = size(xDisp,1);
% apply the margins
ex=ex-sideMargin;
ey=ey-bottomMargin;
%initial index for choosing data on x- and y-axes
sx = 1+sideMargin;
sy = 1+topMargin;
% vector of the selected indexes
ax=sx:reduce:ex;
ay=sy:reduce:ey;

if (elType==606||elType==608)&&(elOrder>=3)% check that there is the proper number of points
  si=size(ax,2);
  a=mod(si-1,elOrder-1);
  ax=ax(1,1:1:si-a);
  ex=ax(end);
  si=size(ay,2);
  a=mod(si-1,elOrder-1);
  ay=ay(1,1:1:si-a);
  ey=ay(end);
  clear a; clear si;
end

% Normalizing measurement with its absolute minimum (useful for linear reconstruction)
%for i=1:1:size(Ln,2)
%  factor = max(max(abs(yDisp(:,:,Ln)),[],1),[],2);    
%  xDisp(:,:,Ln)=xDisp(:,:,Ln)/factor;
%  xDisp(:,:,Ln)=yDisp(:,:,Ln)/factor;
%end


if false% totofalse
  % plot the displacements along some lines in the domain and, if chosen, smooth some measured displacements
  %  with a laplace filter for instance

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
   figure
   subplot(1,2,1)
   surf(xCoor_mod,yCoor_mod,xDisp_mod(:,:,lln),'edgealpha',0);
   view(0,90);
   axis equal;
   title('Initial data: lateral displacement');
   subplot(1,2,2)
   surf(xCoor_mod,yCoor_mod,yDisp_mod(:,:,lln),'edgealpha',0);
   view(0,90);
   axis equal;
   title('axial displacement');
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
    elseif true% laplace smoothing on the boundaries of the domain
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
end% totofalse

% save a copy of the smoothed displacements
if plotBCs  xDispFilteredSmoothed=xDisp; end

% Reducing xDisp and yDisp arrays
xDisp_mod=xDisp(ay,ax,Ln);
yDisp_mod=yDisp(ay,ax,Ln);
m=size(xDisp_mod,1);
n=size(xDisp_mod,2);
if elType==505
  if exist('pDisp');
    pDisp_mod=pDisp(ay,ax,Ln);
  else
    disp('create_input: Warning: Initializing the pressure to 0 as the data does not contain information about the pressure');
    pDisp_mod=zeros(m,n,size(Ln,2));
  end
end
xCoor_mod=xCoor(ay,ax);
yCoor_mod=yCoor(ay,ax);

% node_coor is an assembling of the node number and its x and y
% coordinates

max_node=m*n;
node_coor=zeros(max_node,3);
node_coor(:,1)=[1:1:max_node]';
nodesNum=zeros([n,m]);
nodesNum(:)=node_coor(:,1);

xCoor_mod=xCoor_mod';
yCoor_mod=yCoor_mod';
node_coor(:,2)=xCoor_mod(:);
node_coor(:,3)=yCoor_mod(:);

% Element numbering and its respective nodes
if elType==66||(elType==606&&elOrder==2)||elType==607% making quads
  elem_max=(n-1)*(m-1);
  elem=zeros(elem_max,5);
  elem(:,1)=[1:1:elem_max]';
  tmp=nodesNum(1:1:(n-1),2:1:m);
  elem(:,2)=tmp(:);
  tmp=nodesNum(2:1:n,2:1:m);
  elem(:,3)=tmp(:);
  tmp=nodesNum(1:1:(n-1),1:1:(m-1));
  elem(:,5)=tmp(:);
  tmp=nodesNum(2:1:n,1:1:(m-1));
  elem(:,4)=tmp(:);

  if addElem611
    if ksides==2
      k_poin=zeros(max_node,1);
      kd=1;
      ke=1;
      kc=0;
      bpsides=2*(m-1)-4;
      for k=1:bpsides
        ctr=k+elem_max;
        elem(ctr,1)=ctr;
        %add springs to left,right, and bottom sides (bottom corners ok now!)
        if k<=(m-3)%left side
          elem(ctr,2) = 1+kd*n;
          elem(ctr,3) = 1+(kd+1)*n;
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          kd=kd+1;
        else%right side
          elem(ctr,2) = ke*n+n; 
          elem(ctr,3) = (ke+1)*n+n; 
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          ke=ke+1;
        end
      end
    end
    if ksides==3
      k_poin=zeros(max_node,1);
      kd=1;
      ke=1;
      kc=0;
      bpsides=2*(m-1)-2+1*(n-1);
      for k=1:bpsides
        ctr=k+elem_max;
        elem(ctr,1)=ctr;
        if k<=(m-2)%left side
          elem(ctr,2) = 1+kd*n;
          elem(ctr,3) = 1+(kd+1)*n;
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          kd=kd+1;
        elseif k<=(2*(m-2))%right side
          elem(ctr,2) = ke*n+n; 
          elem(ctr,3) = (ke+1)*n+n; 
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          ke=ke+1;
        else%bottom
          elem(ctr,2) = kc+(m-1)*n+1;
          elem(ctr,3) = elem(ctr,2)+1;
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          kc=kc+1;
        end
      end
    end
    if ksides==4
      k_poin=zeros(max_node,1);
      kd=1;
      ke=1;
      kc=0;
      kf=1;
      bpsides=2*(m-1)-0+2*(n-1);
      for k=1:bpsides
        ctr= k+elem_max;
        elem(ctr,1)=ctr;
        if k<=(m-1)%left side
          elem(ctr,2) = 1+(kd-1)*n;
          elem(ctr,3) = 1+(kd)*n;
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          kd=kd+1;
        elseif k<=(2*(m-1))%right side
          elem(ctr,2) = (ke-1)*n+n; 
          elem(ctr,3) = ke*n+n; 
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          ke=ke+1;
        elseif k<=(2*(m-1)+(n-1))%bottom
          elem(ctr,2) = kc + (m-1)*n + 1;
          elem(ctr,3) = elem(ctr,2) + 1;
          k_poin(elem(ctr,2),1) = 1;
          k_poin(elem(ctr,3),1) = 1;
          kc=kc+1;
        else% top
          elem(ctr,2) = kf;
          elem(ctr,3) = kf+1;
          kf=kf+1;
        end
      end
    end
  end
elseif elType==606&&elOrder==3% make the 9 nodes lagrangian
  elem_max=(n-1)*(m-1)/4;
  elem=zeros(elem_max,10);
  elem(:,1)=[1:1:elem_max]';
  tmp=nodesNum(1:2:(n-2),3:2:m);
  elem(:,2)=tmp(:);
  tmp=nodesNum(3:2:n,3:2:m);
  elem(:,3)=tmp(:);
  tmp=nodesNum(3:2:n,1:2:(m-2));
  elem(:,4)=tmp(:);
  tmp=nodesNum(1:2:(n-2),1:2:(m-2));
  elem(:,5)=tmp(:);
  tmp=nodesNum(2:2:(n-1),3:2:m);
  elem(:,6)=tmp(:);
  tmp=nodesNum(3:2:n,2:2:(m-1));
  elem(:,7)=tmp(:);
  tmp=nodesNum(2:2:(n-1),1:2:(m-2));
  elem(:,8)=tmp(:);
  tmp=nodesNum(1:2:(n-2),2:2:(m-1));
  elem(:,9)=tmp(:);
  tmp=nodesNum(2:2:(n-1),2:2:(m-1));
  elem(:,10)=tmp(:);
elseif elType==608% making quads - different nodes' order
  elem_max=(n-1)*(m-1)/((elOrder-1)^2);
  elem=zeros(elem_max,1+elOrder^2);
  elem(:,1)=[1:1:elem_max]';
  for i=1:1:elOrder
    for j=1:1:elOrder
      tmp=nodesNum(i+[0:(elOrder-1):(n-elOrder)],elOrder+1-j+[0:(elOrder-1):(m-elOrder)]);
      elem(:,1+i+(j-1)*elOrder)=tmp(:);
    end
  end
%  tmp=nodesNum(2:1:n,2:1:m);
%  elem(:,3)=tmp(:);
%  tmp=nodesNum(1:1:(n-1),1:1:(m-1));
%  elem(:,4)=tmp(:);
%  tmp=nodesNum(2:1:n,1:1:(m-1));
%  elem(:,5)=tmp(:);
elseif elType==505% making tris
  disp('elem505, create a mesh that has diagonals in the 2 directions ...for isotropty');
%  return;
  elem_max=2*(n-1)*(m-1);
  elem=zeros(elem_max,4);
  elem(:,1)=[1:1:elem_max]';
  tmp=nodesNum(1:1:(n-1),1:1:(m-1));
  elem(1+2*([1:1:((n-1)*(m-1))]-1),2)=tmp(:);
  tmp=nodesNum(1:1:(n-1),2:1:m);
  elem(1+2*([1:1:((n-1)*(m-1))]-1),4)=tmp(:);
  tmp=nodesNum(2:1:n,1:1:(m-1));
  elem(1+2*([1:1:((n-1)*(m-1))]-1),3)=tmp(:);
  tmp=nodesNum(2:1:n,1:1:(m-1));
  elem(2+2*([1:1:((n-1)*(m-1))]-1),2)=tmp(:);
  tmp=nodesNum(1:1:(n-1),2:1:m);
  elem(2+2*([1:1:((n-1)*(m-1))]-1),4)=tmp(:);
  tmp=nodesNum(2:1:n,2:1:m);
  elem(2+2*([1:1:((n-1)*(m-1))]-1),3)=tmp(:);
else
  disp('elType has a value that is not recognized ... exiting');
  return
end

% Creates a file called input.m (lots of comments are printed for clarity)
if addElem611
  fid = fopen(['kgen' num2str(ksides) '-' filename(1:(end-4)) '.in'], 'wt');
else
  fid = fopen([filename(1:(end-4)) '.in'], 'wt');
end
fprintf(fid,['%% this is an input file for NLACE created on ' datestr(clock) '\n']);
fprintf(fid,['%%  from ' filepath filename '\n']);
fprintf(fid,'\n');
fprintf(fid,'%% NLACE is based on an optimization algorithm (L-BFGS-B)\n');
fprintf(fid,'%%  and a solver (Non-Linear elastic solver)\n');
fprintf(fid,'\n');
if size(Ln,2)==1
  fprintf(fid,['%% In this input file, there is 1 measurement of displacements\n']);
else
  fprintf(fid,['%% In this input file, there are ' num2str(size(Ln,2)) ' different measurements of displacements\n']);
end
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
fprintf(fid,'%% NAME OF THE OUTPUT FILES (PREFIX ONLY - OPTIONAL)\n');
fprintf(fid,'output\n');
if addElem611
  fprintf(fid,['kgen' num2str(ksides) '-' filename(1:(end-4)) '.out\n']);
else
  fprintf(fid,[filename(1:(end-4)) '.out\n']);
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% PRINTING TIMING INFO (FOR PROFILING - OPTIONAL)\n');
fprintf(fid,'%%timingPrint\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% MATERIAL PROPERTIES\n');
if elType==66
  fprintf(fid,'%%  nGaussPtsPerDir regType shearModReg shearRegTVDcst propSetType\n');
  if addElem611
    fprintf(fid,'%%  nGaussPtsPerDir xk11 xk12 xk21 xk22\n');
  end
  fprintf(fid,'materialProperties\n');
  %info3=[3  21  0.0001  6.0d-03 1];
  info3=[3  21  0.00001  1.0d-04 2];
  fprintf(fid,' %d %d %e %e %d\n',info3);
  if addElem611
    info32=[2 1.0d+02 0.000 0.000 5.0d+04];
    fprintf(fid,' %d %e %e %e %e\n' ,info32);
  end
elseif elType==505
  fprintf(fid,'%%  nGaussPtsPerDir regType gammaReg shearModReg gammaTVDcst shearRegTVDcst stabFactor activeStab\n');
  fprintf(fid,'materialProperties\n');
  info3=[3  21  0.0000 0.000000000001  1.0d-04 1.0d-04 0.001 0.001 1.0 1];
  fprintf(fid,' %d %d %e %e %e %e %e %e %f %d\n',info3);
elseif elType==606
  fprintf(fid,'%%  nGaussPtsPerDir regType gammaReg shearModReg gammaTVDcst shearRegTVDcst gammaThres shearModThres\n');
  fprintf(fid,'materialProperties\n');
  info3=[3  21  0.000 0.00001 1.0e-004 1.0d-04 0.01 0.01];
  fprintf(fid,' %d %d %e %e %e %e %e %e\n',info3);
elseif elType==607||elType==608
  fprintf(fid,'%%  nGaussPtsPerDir regType gammaReg shearModReg gammaTVDcst shearRegTVDcst weightL2 weightHx weightHy\n');
  fprintf(fid,'materialProperties\n');
  info3=[3  21  0.000 0.00001 1.0e-004 1.0d-04 1.0 0.0 0.0];
  fprintf(fid,' %d %d %e %e %e %e %e %e %e\n',info3);
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% OPTIONS FOR THE SOLVER (OPTIONAL)\n');
fprintf(fid,'%%  - nCores  : number of cores used by the solver\n');
fprintf(fid,'%%  - solveMethod: 1-> direct solver (default); 2-> iterative solver\n');
fprintf(fid,'%% nCores  solveMethod\n');
fprintf(fid,'soloptions\n');
fprintf(fid,'%d  %d\n',1,1);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% TYPE OF PROBLEM (OPTIONAL)\n');
fprintf(fid,'%% - pbtype: 1-> inverse problem (default); 2-> direct solve\n');
fprintf(fid,'pbtype\n');
fprintf(fid,'%d\n',1);
fprintf(fid,'\n');
fprintf(fid,'\n');
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
fprintf(fid,'%%  - lsteps : number of steps for loading the material (any values if ok)\n');
fprintf(fid,'%%  - ncontin: number of steps for transitioning material properties (any value is ok)\n');
fprintf(fid,'%% nelem  npoin  ndime  mnode  mdofn  nmat  mprops  nmeas  mpoinbc  tol  lsteps  ncontin\n');
fprintf(fid,'femoptions\n');
% modify to handle multiple frames
if elType==66
  if addElem611
    if ksides == 2
      info1 = [(2*(m-1)-4+(m-1)*(n-1)) max_node   2   4   2   2   5  size(Ln,2)  2*(m+n)-4 1.0e-11  1   1];%sides
    elseif ksides == 3
      info1 = [(2*(m-1)-2+1*(n-1)+(m-1)*(n-1)) max_node  2  4  2  2  5  size(Ln,2)  2*(m+n)-4 1.0e-11  1   1];
    elseif ksides == 4
      info1 = [(2*(m-1)-0+2*(n-1)+(m-1)*(n-1)) max_node  2  4  2  2  5  size(Ln,2)  2*(m+n)-4 1.0e-11  1   1];
    end
    fprintf(fid, '%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
  else
    info1 = [(m-1)*(n-1) max_node   2     4     2     1     5  size(Ln,2)  2*(m+n)-4 1.0e-11 1   1];
    fprintf(fid, '%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
  end
elseif elType==505
  info1 = [2*(m-1)*(n-1) max_node   2     3     3     1     10  size(Ln,2)  2*(m+n)-4 1.0e-11 5   1];
  fprintf(fid, '%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
elseif elType==606
  if elOrder==2
    info1 = [(m-1)*(n-1) max_node   2     4     2     1     8  size(Ln,2)  2*(m+n)-4 1.0e-11  5   1];
  elseif elOrder==3
    info1 = [(m-1)*(n-1)/4 max_node   2     9     2     1     8  size(Ln,2)  2*(m+n)-4 1.0e-11  5   1];
  else
    disp('pb with the elOrder...exiting');
    return
  end
  fprintf(fid, '%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
elseif elType==607
  info1 = [(m-1)*(n-1) max_node   2     4     2     1     9  size(Ln,2)  2*(m+n)-4 1.0e-11  5   1];
  fprintf(fid, '%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
elseif elType==608
  info1 = [elem_max max_node   2     elOrder^2     2     1     9  size(Ln,2)  2*(m+n)-4 1.0e-11  5   1];
  fprintf(fid, '%d %d %d %d %d %d %d %d %d %e %d %d\n',info1);
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% BOUNDS FOR THE BOUND-CONSTRAINED OPTIMIZATION (or BOX-CONSTRAINED)\n');
fprintf(fid,'%%  - lpa1: lower bound for parameter 1\n');
fprintf(fid,'%%  - upa1: upper bound for parameter 1\n');
fprintf(fid,'%%  - lpa2: lower bound for parameter 2\n');
fprintf(fid,'%%  - upa2: upper bound for parameter 2, etc\n');
fprintf(fid,'%% lpa1  upa1\n');
fprintf(fid,'%% lpa2  upa2 (there should be mnty lines)\n');
fprintf(fid,'optbounds\n');
%fprintf(fid,'%f %f\n',0.01,4.0);
%fprintf(fid,'%f %f\n',1.0,8.0);
fprintf(fid,'%f %f\n',0.01,10.0);
fprintf(fid,'%f %f\n',1.0,20.0);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% OPTIONS FOR THE BOUND-CONSTRAINED OPTIMIZATION (or BOX-COSNTRAINED)\n');
fprintf(fid,'%%  - iopt   : type of optimization (2->L-BFGS-B; 4->ASA_CG; 5->gencan)\n');
fprintf(fid,'%%  - niter  : maximum number of iterations\n');
fprintf(fid,'%%  - bfgsM  : value of parameter m in L-BFGS-B (usually 10<m<25)\n');
fprintf(fid,'%%  - noutput: interval for saving data\n');
fprintf(fid,'%%  - mnty   : maximum number of fields (types of variables) to be optimized (cf optbounds)\n');
fprintf(fid,'%% iopt  niter  bfgsM  noutput  mnty\n');
fprintf(fid,'optoptions\n');
%fprintf(fid,'%d %d %d %d %d\n',2,20,10,10,2);
fprintf(fid,'%d %d %d %d %d\n',2,20,10,20,2);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% ELEMENT SETS\n');
fprintf(fid,'%% initialElement  finalElement  elementType  numberOfNodesPerElement\n');
fprintf(fid,'elementSet(s)\n');
if elType==66
  info2=[1 (m-1)*(n-1) 66 4];
  fprintf(fid, '%d %d %d %d\n',info2);
elseif elType==505
  info2=[1 2*(m-1)*(n-1) 505 3];
  fprintf(fid, '%d %d %d %d\n',info2);
elseif elType==606
  if elOrder==2
    info2=[1 (m-1)*(n-1) 606 4];
  elseif elOrder==3
    info2=[1 (m-1)*(n-1)/4 606 9];
  end
  fprintf(fid, '%d %d %d %d\n',info2);
elseif elType==607
  info2=[1 (m-1)*(n-1) 607 4];
  fprintf(fid, '%d %d %d %d\n',info2);
elseif elType==608
  info2=[1 elem_max 608 elOrder^2];
  fprintf(fid, '%d %d %d %d\n',info2);
end
if  addElem611 
  e1end=(m-1)*(n-1);
  if ksides==3
    e2end=2*(m-1)-2 + (m-1)*(n-1) + 1*(n-1);
  elseif ksides==2
    e2end=2*(m-1)-4 + (m-1)*(n-1); 
  elseif ksides==4
    e2end=2*(m-1)-0 + (m-1)*(n-1) + 2*(n-1); 
  end
  info23=[(e1end+1) e2end 611 2];
  fprintf(fid, '%d %d %d %d\n',info23);
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% RESTART OPTION (commented by default)\n');
fprintf(fid,'%%restart\n');
fprintf(fid,'uncommentRestartKeywordAndPutYourFileNameHere\n');
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'%% NODES COORDINATES\n');
fprintf(fid,'%% nodeNumber x-coordinate y-coordinate\n');
fprintf(fid,'coordinates\n');
fprintf(fid,'%d %12.6e %12.6e\n',node_coor');
fprintf(fid,'\n');
fprintf(fid,'\n');
% Writing connectivity data
fprintf(fid,'%% ELEMENTS TOPOLOGY\n');
fprintf(fid,'%% elementNumber nodeNumbers ...\n');
fprintf(fid,'connectivity\n');
if size(elem,2)==4
  elem(:,[3,4])=elem(:,[4,3]);% comment or uncomment if need be
  fprintf(fid,' %d   %d   %d   %d\n',elem');
elseif size(elem,2)==5
  if addElem611
    fprintf(fid,' %d   %d   %d   %d   %d\n',elem(1:e1end,1:5)');
    fprintf(fid,' %d   %d   %d\n',elem((e1end+1):e2end,1:3)');
  else
    fprintf(fid,' %d   %d   %d   %d   %d\n',elem');
  end
else
  for i=1:1:size(elem,1)
    for j=1:1:size(elem,2)
      fprintf(fid,' %d ',elem(i,j));
    end
    fprintf(fid,'\n');
  end
end
fprintf(fid,'\n');
fprintf(fid,'\n');
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
boun=2*(m+n)-4; % number of nodes on the boundary of the domain
indexBoundi=[ones(1,n),m*ones(1,n),[2:1:(m-1)],[2:1:(m-1)]];% top, bottom and side boundaries
indexBoundj=[[1:1:n],[1:1:n],ones(1,m-2),n*ones(1,m-2)];
inde=[indexBoundi+(indexBoundj-1)*m]';
nodesNum=nodesNum';
nodesNumB=nodesNum(inde);% nodes number of the nodes on the boundary
inde=[indexBoundj+(indexBoundi-1)*n];
xCoor_bou=xCoor_mod(inde);% x coordinnate of the nodes on the boundary
yCoor_bou=yCoor_mod(inde);% y coordinnate
% provide the square of a distance that is used to find neighbours
dist2ref=((max(yCoor_bou(:),[],1)-min(yCoor_bou(:),[],1)+...
           max(xCoor_bou(:),[],1)-min(xCoor_bou(:),[],1))/40)^2;% square of the tolerance distance (in mm^2)
for lln=1:1:size(Ln,2)
  disp(['create_input: Frame ' num2str(Ln(lln))]);
  fprintf(fid,' %d\n',boun);% number of nodes on the boundary
  inde=[indexBoundi+(indexBoundj-1)*m+(lln-1)*m*n]';
  if addElem611
    disp(['create_input: Imposing spring BCs']);
    %  bound_x_y=[nodesNumB,0*ones(boun,1),xDisp_mod(inde),0*ones(boun,1),yDisp_mod(inde)];
    if ksides == 2
      bound_x_y=[nodesNumB,[1*ones(2*n,1) ; 0*ones((boun - 2*n),1)],xDisp_mod(inde), ...
           [1*ones(2*n,1) ; 0*ones((boun - 2*n),1)],yDisp_mod(inde)];%top & bottom fixed, springs on sides 
    elseif ksides == 3
      bound_x_y=[nodesNumB,[1*ones(n,1) ; 0*ones((boun - n),1)],xDisp_mod(inde), ...
           [1*ones(n,1) ; 0*ones((boun-n),1)],yDisp_mod(inde)];  %top fixed, springs on all other sides 
    elseif ksides == 4
      bound_x_y=[nodesNumB,[0*ones(2*n,1) ; 0*ones((boun - 2*n),1)],xDisp_mod(inde), ...
           [0*ones(2*n,1) ; 0*ones((boun - 2*n),1)],yDisp_mod(inde)];% springs on all 4 sides
    end          
    fprintf(fid, '  %d\n   %d %11.6f  %d %11.6f\n',bound_x_y');
  elseif elType==505 % plane strain
    disp('create_input: plane strain elements, DBC on axial, NBC on lateral and pressure');
    bound_x_y=[nodesNumB,0*ones(boun,1),xDisp_mod(inde),1*ones(boun,1),yDisp_mod(inde),...
               0*ones(boun,1),pDisp_mod(inde)];
    % to avoid rigid body modes, use DBC for a point close to the center of the edge of the boundary closest
    %  to the transducer; average the displacements and find the point closest to its corresponding avg
    % comment out to avoid averaging??
    inds=floor(n/2)+[-floor(n/5):1:ceil(n/5)];
    nptstmp=1;
    indsL=[(inds(1)+[-nptstmp:1:-1]),inds,(inds(end)+[1:1:nptstmp])];% add 3 points on each sides
    avgDsp=zeros(size(inds));% store the average of the lat. disp.
    for i=1:1:size(inds,2)
      avgDisp(1,i)=mean(xDisp_mod(1,indsL(i:1:(i+2*nptstmp)),lln));
    end
      tmp=abs(xDisp_mod(1,inds,lln)-avgDsp);
    inds2=find(tmp==min(tmp(:),[],1));
    inds2=inds2(1);
    bound_x_y(inds(inds2),2)=1;
    fprintf(fid, '  %d\n   %d %11.6f   %d %11.6f   %d %11.6f\n',bound_x_y');
    if plotBCs  boundx{lln}=bound_x_y(:,3); end
  elseif elType==606||elType==66||elType==607||elType==608% plane stress
    % dirichlet BC everywhere for axial disp and homogeneous Neuman everywhere for lateral disp:
    bound_x_y=[nodesNumB,zeros(boun,1),xDisp_mod(inde),zeros(boun,1),yDisp_mod(inde)];
    if false
      disp(['create_input: Lateral disp: interpolating between well-chosen points']);
      % to avoid rigid body modes, fix one or several lateral displacement through dirichlet BC:
      % choose the coordinates of points which are close to a node where the lateral displacement
      %  is enforced with a dirichlet condition
  % note: we could choose the points with the minimal displacements as the homogeneous 
  %  neuman will push the boundaries "outside" the domain...
      if true% fix 17 points with DBC close 
      fixedPoints=[xCoor_mod(floor(n/2),1),yCoor_mod(floor(n/2),1);...
                   xCoor_mod(1,1),yCoor_mod(1,1);...
      	           xCoor_mod(n,1),yCoor_mod(n,1);
		   xCoor_mod(1,floor(m/4)),yCoor_mod(1,floor(m/4));...
		   xCoor_mod(n,floor(m/4)),yCoor_mod(n,floor(m/4));...
		   xCoor_mod(1,floor(3*m/8)),yCoor_mod(1,floor(3*m/8));...
		   xCoor_mod(n,floor(3*m/8)),yCoor_mod(n,floor(3*m/8));...
		   xCoor_mod(1,floor(m/2)),yCoor_mod(1,floor(m/2));...
		   xCoor_mod(n,floor(m/2)),yCoor_mod(n,floor(m/2));...
		   xCoor_mod(1,floor(5*m/8)),yCoor_mod(1,floor(5*m/8));...
		   xCoor_mod(n,floor(5*m/8)),yCoor_mod(n,floor(5*m/8));...
		   xCoor_mod(1,floor(3*m/4)),yCoor_mod(1,floor(3*m/4));...
		   xCoor_mod(n,floor(3*m/4)),yCoor_mod(n,floor(3*m/4));...
		   xCoor_mod(1,floor(7*m/8)),yCoor_mod(1,floor(7*m/8));...
		   xCoor_mod(n,floor(7*m/8)),yCoor_mod(n,floor(7*m/8));...
		   xCoor_mod(1,floor(19*m/20)),yCoor_mod(1,floor(19*m/20));...
		   xCoor_mod(n,floor(19*m/20)),yCoor_mod(n,floor(19*m/20));...
		   ];
      elseif false% 14 points (only on the side boundaries)
      fixedPoints=[xCoor_mod(1,floor(1*m/20)),yCoor_mod(1,floor(1*m/20));...
		   xCoor_mod(n,floor(1*m/20)),yCoor_mod(n,floor(1*m/20));...
		   xCoor_mod(1,floor(m/6)),yCoor_mod(1,floor(m/6));...
		   xCoor_mod(n,floor(m/6)),yCoor_mod(n,floor(m/6));...
		   xCoor_mod(1,floor(m/3)),yCoor_mod(1,floor(m/3));...
		   xCoor_mod(n,floor(m/3)),yCoor_mod(n,floor(m/3));...
		   xCoor_mod(1,floor(m/2)),yCoor_mod(1,floor(m/2));...
		   xCoor_mod(n,floor(m/2)),yCoor_mod(n,floor(m/2));...
		   xCoor_mod(1,floor(2*m/3)),yCoor_mod(1,floor(2*m/3));...
		   xCoor_mod(n,floor(2*m/3)),yCoor_mod(n,floor(2*m/3));...
		   xCoor_mod(1,floor(5*m/6)),yCoor_mod(1,floor(5*m/6));...
		   xCoor_mod(n,floor(5*m/6)),yCoor_mod(n,floor(5*m/6));...
		   xCoor_mod(1,floor(19*m/20)),yCoor_mod(1,floor(19*m/20));...
		   xCoor_mod(n,floor(19*m/20)),yCoor_mod(n,floor(19*m/20));...
		   ];
      elseif false% 12 points (only on the side boundaries)
      fixedPoints=[xCoor_mod(1,floor(1*m/20)),yCoor_mod(1,floor(1*m/20));...
		   xCoor_mod(n,floor(1*m/20)),yCoor_mod(n,floor(1*m/20));...
		   xCoor_mod(1,floor(m/5)),yCoor_mod(1,floor(m/5));...
		   xCoor_mod(n,floor(m/5)),yCoor_mod(n,floor(m/5));...
		   xCoor_mod(1,floor(2*m/5)),yCoor_mod(1,floor(2*m/5));...
		   xCoor_mod(n,floor(2*m/5)),yCoor_mod(n,floor(2*m/5));...
		   xCoor_mod(1,floor(3*m/5)),yCoor_mod(1,floor(3*m/5));...
		   xCoor_mod(n,floor(3*m/5)),yCoor_mod(n,floor(3*m/5));...
		   xCoor_mod(1,floor(4*m/5)),yCoor_mod(1,floor(4*m/5));...
		   xCoor_mod(n,floor(4*m/5)),yCoor_mod(n,floor(4*m/5));...
		   xCoor_mod(1,floor(19*m/20)),yCoor_mod(1,floor(19*m/20));...
		   xCoor_mod(n,floor(19*m/20)),yCoor_mod(n,floor(19*m/20));...
		   ];
      elseif false% 10 points (only on the side boundaries)
      fixedPoints=[xCoor_mod(1,floor(1*m/20)),yCoor_mod(1,floor(1*m/20));...
		   xCoor_mod(n,floor(1*m/20)),yCoor_mod(n,floor(1*m/20));...
		   xCoor_mod(1,floor(m/4)),yCoor_mod(1,floor(m/4));...
		   xCoor_mod(n,floor(m/4)),yCoor_mod(n,floor(m/4));...
		   xCoor_mod(1,floor(2*m/4)),yCoor_mod(1,floor(2*m/4));...
		   xCoor_mod(n,floor(2*m/4)),yCoor_mod(n,floor(2*m/4));...
		   xCoor_mod(1,floor(3*m/4)),yCoor_mod(1,floor(3*m/4));...
		   xCoor_mod(n,floor(3*m/4)),yCoor_mod(n,floor(3*m/4));...
		   xCoor_mod(1,floor(19*m/20)),yCoor_mod(1,floor(19*m/20));...
		   xCoor_mod(n,floor(19*m/20)),yCoor_mod(n,floor(19*m/20));...
		   ];
      elseif false% 8 points (only on the side boundaries)
      fixedPoints=[xCoor_mod(1,floor(1*m/20)),yCoor_mod(1,floor(1*m/20));...
		   xCoor_mod(n,floor(1*m/20)),yCoor_mod(n,floor(1*m/20));...
		   xCoor_mod(1,floor(m/3)),yCoor_mod(1,floor(m/3));...
		   xCoor_mod(n,floor(m/3)),yCoor_mod(n,floor(m/3));...
		   xCoor_mod(1,floor(2*m/3)),yCoor_mod(1,floor(2*m/3));...
		   xCoor_mod(n,floor(2*m/3)),yCoor_mod(n,floor(2*m/3));...
		   xCoor_mod(1,floor(19*m/20)),yCoor_mod(1,floor(19*m/20));...
		   xCoor_mod(n,floor(19*m/20)),yCoor_mod(n,floor(19*m/20));...
		   ];
      elseif false% 6 points (only on the side boundaries)
      fixedPoints=[xCoor_mod(1,floor(1*m/20)),yCoor_mod(1,floor(1*m/20));...
		   xCoor_mod(n,floor(1*m/20)),yCoor_mod(n,floor(1*m/20));...
		   xCoor_mod(1,floor(m/2)),yCoor_mod(1,floor(m/2));...
		   xCoor_mod(n,floor(m/2)),yCoor_mod(n,floor(m/2));...
		   xCoor_mod(1,floor(19*m/20)),yCoor_mod(1,floor(19*m/20));...
		   xCoor_mod(n,floor(19*m/20)),yCoor_mod(n,floor(19*m/20));...
		   ];
      else
        disp('you should choose a number of points to fix with Dirichlet BCs');
      end
      % find the closest node of the boundary to every fixedPoints
      closestNode=zeros(size(fixedPoints,1),1);
      for i=1:1:size(fixedPoints,1)
        dist2=(fixedPoints(i,1)-xCoor_bou).^2+(fixedPoints(i,2)-yCoor_bou).^2;
        ind=find(dist2==min(dist2(:),[],1));
        closestNode(i)=ind(1);
      end
      % find a set of close nodes of the boundary for every closest node
      % for every set, choose the node whose lateral displacement is the closest to the average lateral displacement of the set
      disp(['create_input: Lateral disp: DBC for a subset of ' num2str(size(fixedPoints,1)) ' points on the boundary']);
      for i=1:1:size(closestNode,1)
        dist2=(xCoor_bou(closestNode(i))-xCoor_bou).^2+(yCoor_bou(closestNode(i))-yCoor_bou).^2;
        inds=find(dist2<dist2ref);
        disp_tmp=xDisp_mod(inde(inds));
        avgDisp_tmp=sum(disp_tmp(:),1)/max(size(disp_tmp),[],2);
        disp_tmp=abs(disp_tmp-avgDisp_tmp);
        ind=find(disp_tmp==min(disp_tmp(:),[],1));
        bound_x_y(inds(ind(1)),2)=1;
      end
    elseif false% fix the DBC by hand...
      %inds=floor(n/2)+[-floor(n/5):1:floor(n/5)];
      %tmp=xDisp_mod(1,inds,lln);
      %avgtmp=sum(tmp(:),1)/max(size(tmp),[],2);
      %tmp=abs(tmp-avgtmp);
      %inds2=find(tmp==min(tmp(:),[],1));
      %inds2=inds2(1);
      %bound_x_y(inds(inds2),2)=1;
      %bound_x_y(5,2)=1;
      %bound_x_y(n-5,2)=1;
      %bound_x_y(2*n+floor(m/2),2)=1;
      %bound_x_y(2*n+m+floor(m/2),2)=1;
      %bound_x_y(2*n+floor(3*m/4),2)=1;
      %bound_x_y(2*n+m+floor(3*m/4),2)=1;
    elseif false
      disp('create_input: Lateral disp: Smoothing with linear interpolation between well-chosen points');
      nPtTop=5;% number of points to define segments on the top bc
      nPtBot=5;% number of points to define segments on the bottom bc
      nPtSid=6;% number of points to define segments on the side bc
      % fix a certain number of points and interpolate linearly in-between
      fixedPoints=[xCoor_mod(1,1),yCoor_mod(1,1);...
      	           xCoor_mod(n,1),yCoor_mod(n,1);...
      	           xCoor_mod(1,1),yCoor_mod(1,m);...
      	           xCoor_mod(n,m),yCoor_mod(n,m);...
		   xCoor_mod(floor([1:1:nPtTop]*n/(nPtTop+1)),1),yCoor_mod(floor([1:1:nPtTop]*n/(nPtTop+1)),1);...
		   xCoor_mod(floor([1:1:nPtBot]*n/(nPtBot+1)),m),yCoor_mod(floor([1:1:nPtBot]*n/(nPtBot+1)),m);...
		   xCoor_mod(1,floor([1:1:nPtSid]*m/(nPtSid+1)))',yCoor_mod(1,floor([1:1:nPtSid]*m/(nPtSid+1)))';...
		   xCoor_mod(n,floor([1:1:nPtSid]*m/(nPtSid+1)))',yCoor_mod(n,floor([1:1:nPtSid]*m/(nPtSid+1)))';...
		   ];
      % find the closest node of the boundary to every fixedPoints
      closestNode=zeros(size(fixedPoints,1),1);
      for i=1:1:size(fixedPoints,1)
        dist2=(fixedPoints(i,1)-xCoor_bou).^2+(fixedPoints(i,2)-yCoor_bou).^2;
        ind=find(dist2==min(dist2(:),[],1));
        closestNode(i)=ind(1);
      end
      % find a set of close nodes of the boundary for every closest node
      % for every set, choose the node whose lateral displacement is the closest to the average lateral displacement of the set
      indexAngle=zeros(1,4);
      for i=1:1:4% pour les 4 angles modifier la valeur des angles dans bound_x_y
        dist2=(xCoor_bou(closestNode(i))-xCoor_bou).^2+(yCoor_bou(closestNode(i))-yCoor_bou).^2;
        inds=find(dist2<dist2ref);
        disp_tmp=xDisp_mod(inde(inds));
        avgDisp_tmp=sum(disp_tmp(:),1)/max(size(disp_tmp),[],2);
        disp_tmp=abs(disp_tmp-avgDisp_tmp);
        ind=find(disp_tmp==min(disp_tmp(:),[],1));
        indexAngle(1,i)=inds(ind(1));
      end
      bound_x_y(1,3)=bound_x_y(indexAngle(1,1),3);%bound_x_y(1,5)=bound_x_y(indexAngle(1,1),5);
      bound_x_y(n,3)=bound_x_y(indexAngle(1,2),3);%bound_x_y(n,5)=bound_x_y(indexAngle(1,2),5);
      bound_x_y(n+1,3)=bound_x_y(indexAngle(1,3),3);%bound_x_y(n+1,5)=bound_x_y(indexAngle(1,3),5);
      bound_x_y(2*n,3)=bound_x_y(indexAngle(1,4),3);%bound_x_y(2*n,5)=bound_x_y(indexAngle(1,4),5);
      % top boundary
%disp('top boundary');
      indexAngle=zeros(1,nPtTop+2);
      indexAngle(1,1)=1;
      indexAngle(1,nPtTop+2)=n;
      for i=1:1:nPtTop% pour les points qui ne sont pas aux angles
        j=4+i;% index shift
        dist2=(xCoor_bou(closestNode(j))-xCoor_bou).^2+(yCoor_bou(closestNode(j))-yCoor_bou).^2;
        inds=find(dist2<dist2ref);
        disp_tmp=xDisp_mod(inde(inds));
        avgDisp_tmp=sum(disp_tmp(:),1)/max(size(disp_tmp),[],2);
        disp_tmp=abs(disp_tmp-avgDisp_tmp);
        ind=find(disp_tmp==min(disp_tmp(:),[],1));
        indexAngle(1,i+1)=inds(ind(1));
      end
%indexAngle
      for i=1:1:(nPtTop+1)% number of intervals
	widt=indexAngle(1,i+1)-indexAngle(1,i);
        for j=(indexAngle(1,i)+1):1:(indexAngle(1,i+1)-1)
          bound_x_y(j,3)=((indexAngle(1,i+1)-j)*bound_x_y(indexAngle(1,i),3)+...
	                  (j-indexAngle(1,i))*bound_x_y(indexAngle(1,i+1),3))/widt;
	end
%bound_x_y(indexAngle(1,i):1:indexAngle(1,1+i),3)'
%bound_x_y((indexAngle(1,i)+1):1:indexAngle(1,1+i),3)'-bound_x_y(indexAngle(1,i):1:(indexAngle(1,1+i)-1),3)'
      end
      % bottom boundary
%disp('bottom boundary');
      indexAngle=zeros(1,nPtBot+2);
      indexAngle(1,1)=n+1;
      indexAngle(1,nPtBot+2)=2*n;
      for i=1:1:nPtBot% pour les points qui ne sont pas aux angles
        j=4+nPtTop+i;% index shift
        dist2=(xCoor_bou(closestNode(j))-xCoor_bou).^2+(yCoor_bou(closestNode(j))-yCoor_bou).^2;
        inds=find(dist2<dist2ref);
        disp_tmp=xDisp_mod(inde(inds));
        avgDisp_tmp=sum(disp_tmp(:),1)/max(size(disp_tmp),[],2);
        disp_tmp=abs(disp_tmp-avgDisp_tmp);
        ind=find(disp_tmp==min(disp_tmp(:),[],1));
        indexAngle(1,i+1)=inds(ind(1));
      end
%indexAngle
      for i=1:1:(nPtBot+1)% number of intervals
	widt=indexAngle(1,i+1)-indexAngle(1,i);
        for j=(indexAngle(1,i)+1):1:(indexAngle(1,i+1)-1)
          bound_x_y(j,3)=((indexAngle(1,i+1)-j)*bound_x_y(indexAngle(1,i),3)+...
	                  (j-indexAngle(1,i))*bound_x_y(indexAngle(1,i+1),3))/widt;
	end
%bound_x_y(indexAngle(1,i):1:indexAngle(1,1+i),3)'
%bound_x_y((indexAngle(1,i)+1):1:indexAngle(1,1+i),3)'-bound_x_y(indexAngle(1,i):1:(indexAngle(1,1+i)-1),3)'
      end
      % side boundary 1
%disp('side boundary 1');
      indexAngle=zeros(1,nPtSid+2);
      indexAngle(1,1)=1;
      indexAngle(1,nPtSid+2)=n+1;
      indexAn=zeros(1,nPtSid+2);% for the weight as the index for the side boundaries are not continuous
      indexAn(1,1)=2*n;
      indexAn(1,nPtSid+2)=2*n+m-1;
      for i=1:1:nPtSid% pour les points qui ne sont pas aux angles
        j=4+nPtTop+nPtBot+i;% index shift
        dist2=(xCoor_bou(closestNode(j))-xCoor_bou).^2+(yCoor_bou(closestNode(j))-yCoor_bou).^2;
        inds=find(dist2<dist2ref);
        disp_tmp=xDisp_mod(inde(inds));
        avgDisp_tmp=sum(disp_tmp(:),1)/max(size(disp_tmp),[],2);
        disp_tmp=abs(disp_tmp-avgDisp_tmp);
        ind=find(disp_tmp==min(disp_tmp(:),[],1));
        indexAngle(1,i+1)=inds(ind(1));
	indexAn(1,i+1)=inds(ind(1));
      end
%indexAngle
%indexAn
      for i=1:1:(nPtSid+1)% number of intervals
	widt=indexAn(1,i+1)-indexAn(1,i);
        for j=(indexAn(1,i)+1):1:(indexAn(1,i+1)-1)
          bound_x_y(j,3)=((indexAn(1,i+1)-j)*bound_x_y(indexAngle(1,i),3)+...
	                  (j-indexAn(1,i))*bound_x_y(indexAngle(1,i+1),3))/widt;
	end
%(indexAn(1,i)+1):1:(indexAn(1,i+1)-1)
%bound_x_y(indexAn(1,i):1:indexAn(1,1+i),3)'
%bound_x_y((indexAn(1,i)+1):1:indexAn(1,1+i),3)'-bound_x_y(indexAn(1,i):1:(indexAn(1,1+i)-1),3)'
      end
      % side boundary 2
%disp('side boundary 2');
      indexAngle=zeros(1,nPtSid+2);
      indexAngle(1,1)=n;
      indexAngle(1,nPtSid+2)=2*n;
      indexAn=zeros(1,nPtSid+2);% for the weight as the index for the side boundaries are not continuous
      indexAn(1,1)=2*n+m-2;
      indexAn(1,nPtSid+2)=2*n+2*(m-2)+1;
      for i=1:1:nPtSid% pour les points qui ne sont pas aux angles
        j=4+nPtTop+nPtBot+nPtSid+i;% index shift
        dist2=(xCoor_bou(closestNode(j))-xCoor_bou).^2+(yCoor_bou(closestNode(j))-yCoor_bou).^2;
        inds=find(dist2<dist2ref);
        disp_tmp=xDisp_mod(inde(inds));
        avgDisp_tmp=sum(disp_tmp(:),1)/max(size(disp_tmp),[],2);
        disp_tmp=abs(disp_tmp-avgDisp_tmp);
        ind=find(disp_tmp==min(disp_tmp(:),[],1));
        indexAngle(1,i+1)=inds(ind(1));
	indexAn(1,i+1)=inds(ind(1));
      end
%indexAngle
%indexAn
      for i=1:1:(nPtSid+1)% number of intervals
%indexAngle([i,i+1])
%indexAn([i,i+1])
	widt=indexAn(1,i+1)-indexAn(1,i);
        for j=(indexAn(1,i)+1):1:(indexAn(1,i+1)-1)
          bound_x_y(j,3)=((indexAn(1,i+1)-j)*bound_x_y(indexAngle(1,i),3)+...
	                  (j-indexAn(1,i))*bound_x_y(indexAngle(1,i+1),3))/widt;
	end
%(indexAn(1,i)+1):1:(indexAn(1,i+1)-1)
%bound_x_y(indexAn(1,i):1:indexAn(1,1+i),3)'
      end
      %dbcIndex=1:1:(2*n+2*(m-2));disp('create_input: Lateral disp: DBC everywhere');
      %dbcIndex=(2*n+1):1:(2*n+2*(m-2)); disp('create_input: Lateral disp: DBC on the sides');
      %dbcIndex=(n+1):1:(2*n); disp('create_input: Lateral disp: DBC on the bottom');
      %dbcIndex=zeros(1,0); disp('create_input: Lateral disp: no DBC imposed; will be NBC');
      %dbcIndex=1:1:n; disp('create_input: Lateral disp: DBC on the top');
      dbcIndex=1:1:(2*n); disp('create_input: Lateral disp: DBC on the top and bottom');
      bound_x_y(dbcIndex,2)=1;% Dirichlet for lateral displacements
    else% don't filter the lateral displacements to impose as BC
      %disp('create_input: you should fix at least one Dirichlet BC for lateral disp to avoid rigid body modes');
      disp('create_input: No smoothing of the lateral disp for the bc');
      dbcIndex=zeros(1,0); disp('create_input: Lateral disp: no DBC imposed; will be NBC');
      %dbcIndex=1:1:(2*n+2*(m-2));disp('Lateral disp: DBC everywhere');
      bound_x_y(dbcIndex,2)=1;% Dirichlet for lateral displacements
    end
    if plotBCs  boundx{lln}=bound_x_y(:,3); end
    if sum(bound_x_y(:,2))<0.5
      disp('create_input: Fixing one point on the BC automatically to avoid RBMs on lat. disp.');
      % to avoid rigid body modes, use DBC for a point close to the center of the edge of the boundary closest
      %  to the transducer; average the displacements and find the point closest to its corresponding avg
      inds=floor(n/2)+[-floor(n/5):1:ceil(n/5)];
      indsL=[(inds(1)+[-3:1:-1]),inds,(inds(end)+[1:1:3])];% add 3 points on each sides
      avgDsp=zeros(size(inds));% store the average of the lat. disp.
      for i=1:1:size(inds,2)
        avgDsp(1,i)=mean(xDisp_mod(1,indsL(i:1:(i+6)),lln));
      end
      tmp=abs(xDisp_mod(1,inds,lln)-avgDsp);
      inds2=find(tmp==min(tmp(:),[],1));
      inds2=inds2(1);
      dbcIndex=inds(inds2);
      bound_x_y(dbcIndex,2)=1;% Dirichlet for lateral displacements
    end
    dbcIndex=1:1:(2*n+2*(m-2));disp('create_input: Axial disp: DBC everywhere');
    bound_x_y(dbcIndex,4)=1;% dirichlet for axial displacements
    fprintf(fid, '  %d\n   %d %11.6f   %d %11.6f\n',bound_x_y');
  else
    disp('weird element type ... exiting');
    return
  end
end
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
if addElem611
  datn=4;
else
  datn=2;
end
fprintf(fid,' %d\n',datn);
data(1:max_node,1)=0;% active/inactive optimization parameter 1
if elType==66
  data(1:max_node,2)=2;% initial value parameter 1
else
  data(1:max_node,2)=0.5;%2;% initial value parameter 1
end
data(1:max_node,3)=1;% active/inactive optimization parameter 2
data(1:max_node,4)=1.05;% initial value parameter 2
if addElem611
  xDisp_temp = xDisp_mod(:,:,1)';
  yDisp_temp = yDisp_mod(:,:,1)';
  x_uu = xDisp_temp(:);
  y_uu = yDisp_temp(:);
  for inode = 1:max_node
    data(inode,5)=0*k_poin(inode) ;
    data(inode,7)=0*k_poin(inode) ;%optimization parameter
    data(inode,6) = x_uu(inode);
    data(inode,8) = y_uu(inode);%boundary data (initially == u_meas)
  end
  fprintf(fid,'  %d %11.6f  %d %11.6f  %d %11.6f  %d %11.6f\n',data');
else
  fprintf(fid,'  %d %11.6f  %d %11.6f\n',data');
end


% measured x- and y-displacement & pressure for all the nodes, and weight in the objective function
real_disp=zeros(max_node,6);
if false% weight of the lateral displacement in the objective function
  % make a variable map depending on the position of the point
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
real_disp(:,3)=1.0;
real_disp(:,5)=0.0;

if false% variable weight in the objective function depending on the compression level
  avgAxialDisp=zeros(size(Ln,2),1);
  if false% use the average value of the max compression
    disp('create_input: Weight in the obj func equals the max compression');
    for lln=1:1:size(Ln,2)
      avgAxialDisp(lln,1)=mean(yDisp_mod(size(yDisp_mod,1),:,lln));
    end
    maxAvgAxialDisp=max(avgAxialDisp,[],1);
    avgAxialDisp=avgAxialDisp/maxAvgAxialDisp;
    for lln=size(Ln,2):-1:2
      avgAxialDisp(lln,1)=avgAxialDisp(lln,1)/avgAxialDisp(lln-1,1);
    end
  elseif false% use the average value of the top and bottom line
    disp('create_input: Weight in the obj func equals the mean of max and min compressions');
    for lln=1:1:size(Ln,2)
      avgAxialDisp(lln,1)=mean([yDisp_mod(size(yDisp_mod,1),:,lln),yDisp_mod(1,:,lln)]);
    end
    maxAvgAxialDisp=max(avgAxialDisp,[],1);
    avgAxialDisp=avgAxialDisp/maxAvgAxialDisp;
    for lln=size(Ln,2):-1:2
      avgAxialDisp(lln,1)=avgAxialDisp(lln,1)/avgAxialDisp(lln-1,1);
    end
  elseif false% use another value for the weights (maxCompression-minCompression)
    for lln=1:1:size(Ln,2)
      avgAxialDisp(lln,1)=mean(yDisp_mod(size(yDisp_mod,1),:,lln))-mean(yDisp_mod(1,:,lln));
    end
    maxAvgAxialDisp=max(avgAxialDisp,[],1);
    avgAxialDisp=avgAxialDisp/maxAvgAxialDisp;
    for lln=size(Ln,2):-1:2
      avgAxialDisp(lln,1)=avgAxialDisp(lln,1)/avgAxialDisp(lln-1,1);
    end
  else % assumes that the noise in the measurement is normal distribution for every incremental frame
    disp('create_input: Weight in obj. func. assuming white Gaussian noise and uniform compression increments');
    disp('                weight by the sqrt of the max compression');
    avgAxialDisp=sqrt(Ln');
    maxAvgAxialDisp=max(avgAxialDisp,[],1);
    avgAxialDisp=avgAxialDisp/maxAvgAxialDisp;
    for lln=size(Ln,2):-1:2
      avgAxialDisp(lln,1)=avgAxialDisp(lln,1)/avgAxialDisp(lln-1,1);
    end
  %else % use sqrt(max-min axial disp.)? see what is more reasonnable
  end
else
  disp('create_input: Every frame is given the same weight in the objective function');
  avgAxialDisp=ones(size(Ln,2),1);% no multiplicative factor
end

fprintf(fid,'\n');
fprintf(fid,'\n');
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
fprintf(fid,'%%  pair is the displacement along the axial direction and, when there is a third pair,\n');
fprintf(fid,'%%  it corresponds to a measure of the pressure inside the material.\n');
fprintf(fid,'measuredDisplacements\n');
for lln=1:1:size(Ln,2)
  uu=xDisp_mod(:,:,lln)';
  vv=yDisp_mod(:,:,lln)';
  real_disp(:,3)=real_disp(:,3)/avgAxialDisp(lln,1);
  real_disp(:,2)=uu(:);
  real_disp(:,4)=vv(:);
  if elType==505
    uu=pDisp_mod(:,:,lln);
    real_disp(:,6)=uu(:);
    fprintf(fid,'%8.3f %11.6f %8.3f %11.6f %8.3f %11.6f\n', real_disp');
  elseif elType==606||elType==66||elType==607||elType==608
    real_disp(:,6)=0.0;
    fprintf(fid,'%8.3f %11.6f %8.3f %11.6f\n', real_disp(:,[1:4])');
  else
    disp('error: weird type of element ... returning!');
    return
  end
end
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'end\n');


fclose(fid);

% if plotBCs
%   xDispi=xDispi(ay,ax,Ln);
%   xDispFiltered=xDispFiltered(ay,ax,Ln);
%   xDispFilteredSmoothed=xDispFilteredSmoothed(ay,ax,Ln);
%   for lln=1:1:size(Ln,2)
%     indexes=indexBoundi+(indexBoundj-1)*m+n*m*(lln-1);
%     boundi{lln}=xDispi(indexes);
%     boundNlm{lln}=xDispFiltered(indexes);
%     boundNlms{lln}=xDispFilteredSmoothed(indexes);
%   end
%   inds=[[1:1:n],[(2*n+m-2+1):1:(2*n+2*(m-2))],[(2*n):-1:(n+1)],[(2*n+m-2):-1:(2*n+1)]];
%   for lln=1:1:size(Ln,2)
%     figure
%     set(gcf,'color','w');
%     hold on
%     plot(boundi{lln}(inds),'k-','linewidth',2);
%     plot(boundNlm{lln}(inds),'b-','linewidth',2);
%     plot(boundNlms{lln}(inds),'r-','linewidth',2);
%     plot(boundx{lln}(inds),'g-','linewidth',2);
%     legend('initial','filtered','filtered+smoothed','imposed');
%     title(['lateral displacements at the boundaries for frame ' num2str(Ln(1,lln))]);
%   end
% end

