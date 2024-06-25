clear all
codes = {'a' 'b' 'c'}; % Кодировка
% n1=1; n2=0; % Счётчики количества буквенных итераций n1-триплетами n2-1:3
n=1;
numtasks=2; %потом заменить 3 на length(tasksA)
% cursor=numtasks;
counters(1:numtasks,1)='a';
counters(numtasks)=counters(numtasks)-1; % Для правильной работы реле в цикле

% for i=1:length(codes):length(codes)^numtasks
while ~prod(counters(1:end)=='c')
    %     if mod(n1+n2, length(codes))==0
    %         N=(n1+n2)/length(codes);
    %         for i=length(codes):length(codes):length(codes)^numtasks
    % %             buf=i;
    %             if buf < gcd(n1+n2, i)
    %                 N=(n1+n2)/
    %                 buf=gcd(n1+n2, i);
    %             end
    %         end
    %     end
    %     if n1>0
    %         for i=1:length(codes)-1
    %             if mod((n1+n2),length(codes)^i)==0
    %                 counters(cursor)=counters(cursor)+1;
    %             end
    %         end
    %     end
    %     n2=n2+1;
    %     if n2==1 %n1mod3==0
    %         counters(end)='a'; % Перещёлкивает реле задач на букву вперёд поочереди справа на лево
    %         cursor_buf=cursor;
    %
    %     else
    %         counters(end)=counters(end)+1;
    %     end

    %         counters(cursor+(n1/length(codes)))=counters(cursor+(n1/length(codes)))+1; % Перещёлкивает реле задач на букву вперёд поочереди справа на лево
    %         counters(n_count+mod(n1,length(codes)))=counters(n_count+mod(n1,length(codes)))+1;
    %     N=divisors(n); % Получаю все делители n
    N=gcd(n,length(codes)^numtasks); % Получаю наибольший общий делитель между кол-вом кодировок и кол-вом задач с альтернативами N==length(codes)^

    cursor=round(log(N)/log(length(codes))); % Позиция в реле для перещёлкивания (обозвать cursor?) / СТЕПЕНЬ КОЛИЧЕСТВА ЭЛЕМЕНТОВ КОДИРОВКИ

    if gcd(n+1,length(codes)^numtasks)<N % Если на этом ходу переключение старшего разряда, а в следующем снова наименьшего...
        cursor_buf=cursor;
        cursor=0;
        counters(length(codes)-cursor)=counters(length(codes)-cursor)+1; % Определение позиции в которой необходима инкрементация реле и перещёлкивание
        flag=true; 
    end %... то меняем пока что младший разряд, но поднимаем флаг для следующей итерации, что необходимо повысить старший разряд

    if ~flag % Штатный режим
        counters(length(codes)-cursor)=counters(length(codes)-cursor)+1; % Определение позиции в которой необходима инкрементация реле и перещёлкивание
    end

    
    if prod(counters(1:end)=='c')
        break
    end

    if flag %&& counters(numtasks-cursor+1)==codes{end} %N>1
        counters(length(codes)-cursor_buf)=counters(length(codes)-cursor_buf)+1;
        
        for i=numtasks-cursor_buf+1:numtasks
            if i==numtasks
            counters(i)='a';
            counters(i)=counters(i)-1;
            else
                counters(i)='a';
            end
        end
        clear cursor_buf
        flag=false;
    end
    %     N=n/



    %     if mod((n1+n2), increm)==0
    %         n2==length(codes)^n1+n2/numtasks
    %
    %         cursor=cursor-1; % Переход на следующие задачи
    %         increm=increm*length(codes);
    %     end
    %         n1=n1+length(codes); n2=0;
    %
    %
    n=n+1;
end

aaa=1