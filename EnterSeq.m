function newseq = EnterSeq(seq,NOC)
lenseq = length(seq);
newlen = 0;
for i = 1:1:lenseq
    switch seq(i)
        case '2'
            lenseq = lenseq + 1;
            newlen = newlen + 1;
        case '1'
            lenseq = lenseq + 1;
            newlen = newlen + 1;
        case '0'
            lenseq = lenseq + 1;
            newlen = newlen + 1;
        case '-'
            lenseq = lenseq + 0;
    end
end
lenseq = lenseq - 1;
k = 1;
minus = 0;
j = 1;
newseq = zeros(1,newlen);
while k <= lenseq
    switch seq(k)
        case '2'
            seq = [seq(1:k),',',seq((k+1):length(seq))];
            k = k + 2;
            if minus == 1
                %disp('-2');
                newseq(j) = -2;
                j = j + 1;
                minus = 0;
            else
                %disp('2');
                newseq(j) = 2;
                j = j + 1;
            end   
        case '1'
            seq = [seq(1:k),',',seq((k+1):length(seq))];
            k = k + 2;
            if minus == 1
                %disp('-1');
                newseq(j) = -1;
                j = j + 1;
                minus = 0;
            else
                %disp('1');
                newseq(j) = 1;
                j = j + 1;
            end        
        case '0'
            seq = [seq(1:k),',',seq((k+1):length(seq))];
            %disp('0');
            newseq(j) = 0;
            j = j + 1;
            k = k + 2;
        case '-'
            minus = 1;
            k = k + 1;
    end
end
newseq = repmat(newseq,1,NOC);
end