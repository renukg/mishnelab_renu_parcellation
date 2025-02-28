function Bou = getBoundariesFromZerosOnes(mat)
NROWS = 32;
NCOLS = 81;
Bou=false([NROWS, NCOLS]);

for k = 1:size(mat,2)
    %[B,~] = bwboundaries(reshape(mat(:,k),[256 256]));
    [B,~] = bwboundaries(reshape(mat(:,k),[32 81]));
    maxsel = -inf;seli = 0;
    for l = 1:length(B)
        if size(B{l},1) > maxsel
            seli = l;
            maxsel=size(B{l},1);
        end
    end
    for i=1:size( B{seli},1)
        Bou(B{seli}(i,1), B{seli}(i,2))=1;
    end
end
Bou = bwskel(Bou); % to make boundaries thinner