
import csv

csv_crispr = file('/home/nonojh/crisprCas9/SSL/data/Crispr_Data.csv','rb')
reader_crispr = csv.reader(csv_crispr)
wcsvfile = file('/home/nonojh/crisprCas9/SSL/data/Crispr_SSL_pre.txt','wb')
writer = csv.writer(wcsvfile)
raw_num = 0
for row_crispr in reader_crispr:
    if row_crispr[0] == 'index':
        continue
    data = [row_crispr[1],row_crispr[2],row_crispr[4],row_crispr[5],row_crispr[6],row_crispr[7],row_crispr[8],
            row_crispr[9],row_crispr[10],row_crispr[11],row_crispr[12],row_crispr[13],row_crispr[14]]
    writer.writerow(data)
    raw_num = raw_num+1;

csv_crispr.close()
wcsvfile.close()
