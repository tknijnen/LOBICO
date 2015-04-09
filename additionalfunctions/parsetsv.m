function [Samples,Features,IC50s,MutationMatrix] = parsetsv(file);

%% check number of tabs
fid = fopen(file);
NL = 10;
tline = cell(NL,1);
not = zeros(NL,1);
for n = 1:10;
    tline{n} = fgetl(fid);
    display(tline{n});
    not(n) = sum(double(tline{n})==9);
end
fclose(fid);

%% Create format string
unot = unique(not);
if length(unot)~=1;
    error(num2str(unot));
else
    rsh = [];
    rsd = [];
    for n = 1:unot+1;
        if n>1;
            rsd = cat(2,rsd,'%n');
        else
            rsd = cat(2,rsd,'%s');
        end
        rsh = cat(2,rsh,'%s');
    end
end

%% read in header
fid = fopen(file);
H = textscan(fid,rsh,1,'Delimiter','\t');
fclose(fid);

%% read in data
fid = fopen(file);
D = textscan(fid,rsd,'Delimiter','\t','Headerlines',1);
fclose(fid);

%% Parse
N = length(D{1});
P = length(H)-2;

Samples = D{1};
IC50s = D{2};

MutationMatrix = NaN(N,P);
Features = cell(1,P);

for p = 1:P;
    Features{p} = H{p+2}{1};
    MutationMatrix(:,p) = D{p+2};
end


    
    
    

