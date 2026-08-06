% This function is used to detect the number of decimal places in number
% "K". For example, D(0.002) = 3, D(3.89) = 2. The number "K" should not
% bigger than 10. The detection limit is 14 decimal places.

function n = decimal_num(K)
formatSpec = '%.14f';
str = sprintf(formatSpec, K);

dotPos = find(str == '.');
if isempty(dotPos)
    n = 0;
    return;
end

decimalPart = str(dotPos+1:end);
decimalPart = regexprep(decimalPart, '0+$', '');
n = length(decimalPart);