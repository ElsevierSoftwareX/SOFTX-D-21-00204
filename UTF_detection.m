function [encoding, bytes_per_char, BOM_size, byte2char] = UTF_detection(FILENAME)
    
    %Detection process
    %if the input file does not exist, no point doing anything else
    
    [fid, message] = fopen(FILENAME, 'r');   %not 'rt' !!
    if fid < 0
        fprintf(2, 'Failed to open file "%s" because: "%s"\n', FILENAME, message);
        encoding = ''; bytes_per_char = 0; BOM_size = 0; byte2char = @(B) '';
        return
    end
    
    %if there are no byte marks, return default encoding
    
    [~, ~, ~, default_encoding] = fopen(fid);
    
    %The information about the encoding of the file can be deduced from the
    %first 4 bytes.
    
    firstbytes = fread(fid, [1,4], '*uint8');
    fclose(fid);
       
    trust_unicode = ~verLessThan('matlab', '7.5');
       
    %in the below, swapbytes is used as needed to change the input order to
    %the order of the current system.
    
    [~, ~, endian] = computer();
    islittle = strcmpi(endian, 'L');
    
    u8pad = @(B, N) reshape([uint8(B(:)); zeros(N-(mod(length(B)-1,N)+1),1,'uint8')],1,[]);
    
    %Encoding process
    
    if length(firstbytes) >= 2 && all(firstbytes(1:2) == [254, 255])   %UTF16BE
        encoding = 'UTF16BE';
        bytes_per_char = 2;
        BOM_size = 2;
        if trust_unicode
            byte2char = @(B) native2unicode(u8pad(B,bytes_per_char), encoding);
        elseif islittle
            byte2char = @(B) char(swapbytes(typecast(u8pad(B,bytes_per_char), 'uint16')));
        else
            byte2char = @(B) char(typecast(u8pad(B,bytes_per_char), 'uint16'));
        end
    elseif length(firstbytes) >= 4 && all(firstbytes(1:4) == [255, 254, 0, 0])    %UTF32LE
        encoding = 'UTF32LE';
        bytes_per_char = 4;
        BOM_size = 4;
        if trust_unicode
            byte2char = @(B) native2unicode(u8pad(B,bytes_per_char), encoding);
        elseif islittle
            byte2char = @(B) char(typecast(u8pad(B,bytes_per_char),'uint32'));
        else
            byte2char = @(B) char(swapbytes(typecast(u8pad(B,bytes_per_char), 'uint32')));
        end
        % warning('32 bit Unicode detected (UTF32LE); MATLAB only stores 16 bits per character.');
    elseif length(firstbytes) >= 2 && all(firstbytes(1:2) == [255, 254])   %UTF16LE
        encoding = 'UTF16LE';
        bytes_per_char = 2;
        BOM_size = 2;
        if trust_unicode
            byte2char = @(B) native2unicode(u8pad(B,bytes_per_char), encoding);
        elseif islittle
            byte2char = @(B) char(typecast(u8pad(B,bytes_per_char),'uint16'));
        else
            byte2char = @(B) char(swapbytes(typecast(u8pad(B,bytes_per_char), 'uint16')));
        end
    elseif length(firstbytes) >= 4 && all(firstbytes(1:4) == [0, 0, 254, 255])   %UTF32BE
        encoding = 'UTF32BE';
        bytes_per_char = 4;
        BOM_size = 4;
        if trust_unicode
            byte2char = @(B) native2unicode(u8pad(B,bytes_per_char), encoding);
        elseif islittle
            byte2char = @(B) char(swapbytes(typecast(u8pad(B,bytes_per_char), 'uint32')));
        else
            byte2char = @(B) char(typecast(u8pad(B,bytes_per_char), 'uint32'));
        end
        % warning('32 bit Unicode detected (UTF32BE); MATLAB only stores 16 bits per character');
    elseif length(firstbytes) >= 3 && all(firstbytes(1:3) == [239, 187, 191])    %UTF8
        encoding = 'UTF8';
        bytes_per_char = 1;
        BOM_size = 3;
        if trust_unicode
            byte2char = @(B) native2unicode(u8pad(B,bytes_per_char), encoding);
        else
            byte2char = @(B) char(u8pad(B,bytes_per_char));
        end
    else
        encoding = default_encoding;
        bytes_per_char = 1;
        BOM_size = 0;
        if trust_unicode
            byte2char = @(B) native2unicode(u8pad(B,bytes_per_char), encoding);
        else
            byte2char = @(B) char(u8pad(B,bytes_per_char));
        end
    end
end