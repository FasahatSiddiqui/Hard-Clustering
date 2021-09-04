%%
%{
file_path='E:\Matlab_code\Segmentation\Images\';
file_name='House.bmp';
file = strcat(file_path,file_name);
cluster_number=5;
Terminate_condition=0.1;
%}
%%
function[center, output_data, labels]=OKM(Terminate_condition, cluster_number, file)
input_data=imread(file);
input_data=rgb2gray(input_data);
tic;
[m,n]=size(input_data);
input_data=double(input_data);
ran=randperm(255);
input_data=input_data(:);
res_data=input_data';
length_data=m*n;
output_data=zeros(1,length_data);
count_it=2;
MSE(1)=0;
for l=1:cluster_number
    center(l,1)=ran(l);  
end
center=sort(center,'ascend');
pixel_dist=zeros(cluster_number,length_data);
while 1
    pts_final=zeros(cluster_number,length_data);
    cluster_ptsdist=zeros(cluster_number,length_data);
    for i=1:cluster_number
        pixel_dist(i,:)=dist(input_data,center(i));     % Distance Cal.
    end
    update_dist=pixel_dist;
    [dist_val,labels]=min(pixel_dist,[],1);
    for i=1:cluster_number
        [r,c]=find(labels==i);
        [b(i),a(i)]=size(c);                            % No. Elements in a region
        if a(i)~=0
            for j=1:a(i)
                cluster_ptsdist(i,j)=dist_val(1,c(j));  % Diff. b/w pixel & its region center
            end
        end
        if a(i)==0
            cluster_ptsdist(i,:)=0;
        end
    end
    fitness=sum(cluster_ptsdist');                % Fitness of cluster
    Q=cluster_ptsdist.^2;
    FVA=sum(Q');
    MSE(count_it)=(sum(FVA(:)))/length_data;      % MSE (mean square error)
    [fit_val,fit_labels]=sort(fitness,'descend'); % Cluster Fitness in descending order
    [E_cluster,E_labels]=find(a==0);              % Find a Null Cluster
    empty_cluster=length(E_labels);
    update_cluster=cluster_number-empty_cluster;
    if empty_cluster>0                            % 1st Condition for dead unit reducing
        Count_em=0;
        center_update=0;
        for i=1:cluster_number
            if a(i)>0
                Count_em=Count_em+1;
                center_update(Count_em,1)=center(i);                % Ci in assending order after drop dead Clu
            end
        end
        Mid_pts=center_update(1:update_cluster-1,1)+center_update(2:update_cluster,1);    % Pts finding 
        Mid_pts=(Mid_pts./2)';                                      % Pts finding 
        for i=1:(update_cluster-1)
            pts_dist=dist(center(fit_labels(i)),Mid_pts);           % dist. b/w pts and high cluster center
            [VDI,IDI]=sort(pts_dist);                               % indices values
            pts(i,1)=Mid_pts(IDI(1));                               % Final available pts
            Mid_pts(IDI(1))=500;
        end 
        length_pts=length(pts);
        if empty_cluster<=length_pts
            for i=1:empty_cluster
                center_update(update_cluster+i,1)=pts(i,1);
            end
        end
        if empty_cluster>length_pts
            for i=1:length_pts
                center_update(update_cluster+i,1)=pts(i,1);         % Pts int. in cluster center
            end
        end
        center_update=sort(center_update);
        center=center_update;
        for i=1:cluster_number
            pixel_dist(i,:)=dist(input_data,center(i));             % Distance Cal.
        end
        update_dist=pixel_dist;
    end
    if empty_cluster==0                                             % 2nd Condition for reduce MSE
        for i=1:(cluster_number-1)
            pts_f(i)=(center(i)+center(i+1))./2;                    % Pts finding 
        end
        for i=1:(cluster_number-1)
            if fitness(i)>fitness(i+1)
                TDD=update_dist(i,:);
                TDD(res_data==pts_f(i))=255;                        % Pts transferring
                update_dist(i,:)=TDD(1,:);
            end
            if fitness(i+1)>fitness(i)
                TDD=update_dist(i+1,:);
                TDD(res_data==pts_f(i))=255;                        % Pts transferring
                update_dist(i+1,:)=TDD(1,:);
            end
        end
    end
    [dist_val,labels]=min(update_dist,[],1);
    for i=1:cluster_number
        [r,c]=find(labels==i);
        [b(i),a(i)]=size(c);                                        % No. of pixel point in group
        if a(i)~=0
            for j=1:a(i)
                pts_final(i,j)=res_data(1,c(j));
            end
            Sum(i,1)=sum(pts_final(i,:));
            mean(i,1)=(Sum(i,1)/a(i));
        end
        if a(i)==0;
            mean(i,1)=0; 
            pts_final(i,:)=0;
        end
    end
    ME=mean;
    check_terminate=((MSE(count_it)-MSE(count_it-1)).^2).^0.5;      % Condition Defining.
    count_it=count_it+1
    if check_terminate < (1/(2*cluster_number)) || check_terminate < Terminate_condition
        break;
    end
    center=mean;
end
Time=toc
center=round(center);
center=center';
output_data=center(labels);
output_data=reshape(output_data,m,n);
output_data=uint8(output_data);
imshow(output_data);
end