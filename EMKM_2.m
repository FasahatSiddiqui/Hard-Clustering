%{
% apply the commented code to set input values for the function
file_path='E:\Matlab_code\Segmentation\Images\';
file_name='House.bmp';
file = strcat(file_path,file_name);
input_data=imread(file);
cluster_number=4;
%}
%%
function[labels,output_data]=EMKM_2(cluster_number, input_data)
count=2;
MSE(1)=0;
input_data=rgb2gray(input_data);
tic;
[m,n]=size(input_data);
KL=0; KT=0;
cluster_number=4;
lenght_data=n*m;
T=1/2; Ta=T; Tb=T; %Terminate criteria
dimage=double(input_data);
pixel_dist=zeros(cluster_number,lenght_data);
ran=randperm(255);
for l=1:cluster_number
    center(l,1)=ran(l);  
end
updated_center=center;
lenght_data=n*m;
dtype_inputdata=dimage(:); res_inputdata=dtype_inputdata';
pts=zeros(cluster_number,lenght_data);
for i=1:cluster_number
    pixel_dist(i,:)=dist(dtype_inputdata,center(i));
end
[S,labels]=min(pixel_dist,[],1);
for i=1:cluster_number
    [r,c]=find(labels==i); 
    [b(i),a(i)]=size(c); 
    if a(i)~=0
        for j=1:a(i)
            pts(i,j)=res_inputdata(1,c(j)); % no. of pixel in each cluster
        end
        Sum(i,1)=sum(pts(i,:));
        mean(i,1)=(Sum(i,1)/a(i));
    end
    if a(i)==0
        pts(i,:)=0;
        mean(i,1)=0; Sum(i,1)=0;
    end
end
updated_center=mean;
N=a;
%%
while 1              %loop for Tb
    NNAD=1; %intensity value add during the data pts transferring process
    while 1          %loop for Ta
        cluster_pts=zeros(cluster_number,lenght_data);
        for i=1:cluster_number
            update_dist(i,:)=dist(dtype_inputdata,updated_center(i));
        end
        [SSK,IIK]=min(update_dist,[],1);
        for i=1:cluster_number
            [rr,cc]=find(IIK==i); 
            [bb(i),aa(i)]=size(cc); % no. of pixel point in group.
            if aa(i)~=0
                for j=1:aa(i)
                    cluster_pts(i,j)=SSK(1,cc(j));
                end
            end
            if aa(i)==0
                cluster_pts(i,:)=0;
            end
        end
        for i=1:cluster_number
            F=sum(cluster_pts(i,:));
            fitness(i,1)=F;
        end
        pt=zeros(cluster_number,lenght_data);
        [VA, IF]=sort(fitness,'descend');
        B_dist=dist(updated_center,updated_center(IF(1)));
        [VAM, IM]=sort(B_dist,'ascend');
        if VA(cluster_number)>=Ta*VA(1)
            True=1    % indicate inside (Ta) loop break
            break;
        end
        NEG=updated_center(IM(2));
        NCC=((NEG-updated_center(IF(1))).^2).^0.5;
        NCF=(NCC/2);
        NCS=(NCC/2)+NNAD;
        [Rr,Cc]=find(update_dist(IF(1),:)<=NCF & update_dist(IM(2),:)<=NCS);
        if length(Cc)~=0
            [Frs,Fcs]=size(Cc);
            for ppt=1:Fcs
                update_dist(IM(2),Cc(ppt))=0;
            end
        end
        NNAC=NCF/2;
        if fitness(IM(2))~=VA(cluster_number)
            B_Dist=dist(updated_center,updated_center(IF(cluster_number)));
            [VAMS, IMS]=sort(B_Dist,'ascend');
            NEGS=updated_center(IMS(2));
            NCCS=((NEGS-updated_center(IF(cluster_number))).^2).^0.5;
            NCFS=(NCCS/2);
            NCSS=(NCCS/2)+NNAD;
            [Rsr,Csr]=find(update_dist(IF(cluster_number),:)<=NCSS & update_dist(IMS(2),:)<=NCFS);
            if length(Csr)~=0
                [Frss,Fcss]=size(Csr);
                for ppt=1:Fcss
                    update_dist(IF(cluster_number),Csr(ppt))=0;
                end
            end
            NNAC=NCFS/2;
        end
        Xi=dtype_inputdata';
        [Sks,Is]=min(update_dist,[],1);
        for i=1:cluster_number
            [RR,CC]=find(Is==i);
            Nc(i)=length(CC);
            for ppt=1:Nc(i)
                pt(i,ppt)=Xi(1,CC(ppt));
            end
        end
        NF=sum(pt');
        NI=NF./Nc;
        NNi=round(NI);
        NNi=NNi';
        NNAD=NNAD+1;
        if NNi ~= updated_center
            NNAD=1;
        end
        if NNAD==NNAC
            NNAD=1;
        end
        updated_center=NNi;
        Ta=Ta-(Ta/cluster_number);
        KL=KL+1 % indicate inside loop complete
    end
    center=updated_center;
    pts=zeros(cluster_number,lenght_data);
    for i=1:cluster_number
        pixel_dist(i,:)=dist(dtype_inputdata,center(i));
    end
    [S,labels]=min(pixel_dist,[],1);
    for i=1:cluster_number
        [r,c]=find(labels==i);
        [b(i),a(i)]=size(c); % no. of pixel point in group.
        if a(i)~=0
            for j=1:a(i)
                pts(i,j)=res_inputdata(1,c(j));
            end
            Sum(i,1)=sum(pts(i,:));
            mean(i,1)=(Sum(i,1)/a(i));
        end
        if a(i)==0
            pts(i,:)=0;
            mean(i,1)=0; Sum(i,1)=0;
        end
    end
    updated_center=mean; 
    N=a;
    if VA(cluster_number)>=Tb*VA(1)
        True=2  % indicate Top loop (Tb) break
        break;        
    end
    KL=0;
    Ta=T;
    Tb=Tb-(Tb/cluster_number);
    KT=KT+1     % indicate Top loop complete
end
Time=toc
updated_center=round(updated_center);
output_data=updated_center(labels);
output_data=reshape(output_data,m,n);
output_data=uint8(output_data);
imshow(output_data)
end