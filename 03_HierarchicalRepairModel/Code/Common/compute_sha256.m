function signature=compute_sha256(data)
if ischar(data) || isstring(data)
    data=unicode2native(char(data),'UTF-8');
end
try
    digest=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    digest.update(typecast(uint8(data(:)),'int8'));
    bytes=typecast(digest.digest(),'uint8');
    signature=lower(reshape(dec2hex(bytes,2).',1,[]));
catch ME
    error('RepairTask:SHA256Unavailable', ...
        'SHA-256 is required for reproducibility: %s',ME.message);
end
end
