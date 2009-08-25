function write_binary( m, filename );
  % write matrix or vector to file in binary format (double)
  fid = fopen(filename, 'wb');
  nel = size(m(:),1);
  fwrite(fid, m(:), 'double');
  fclose(fid);