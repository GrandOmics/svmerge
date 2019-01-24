import sys
import sv_vcf

def remove_chr_in_header(header_lines):
    b37_contigs = {'1':249250621, '2':243199373, '3':198022430, '4':191154276, 
        '5':180915260, '6':171115067, '7':159138663, '8':146364022, '9':141213431,
        '10':135534747, '11':135006516, '12':133851895, '13':115169878, '14':107349540,
        '15':102531392, '16':90354753, '17':81195210, '18':78077248, '19':59128983,
        '20':63025520, '21':48129895, '22':51304566, 'X':155270560, 'Y':59373566, 'MT':16569}
    chr_map = dict([("chr"+str(i), str(i)) for i in range(1,23)])
    chr_map["chrX"] = "X"
    chr_map["chrY"] = "Y"
    # chr_map["chrM"] = "MT"

    header_new = []

    for line in header_lines:
        line = line.strip()
        if line.startswith("##contig="):
            fields = line[10:-1].split(",")
            id_map = dict([(i.split("=")[0].strip(), i.split("=")[1].strip()) for i in fields])
            if id_map["ID"] in chr_map:
                if int(id_map["length"]) == b37_contigs[chr_map[id_map["ID"]]]:
                    id_map["ID"] = chr_map[id_map["ID"]]
                    header_new.append("##contig=<ID={},length={}>".format(id_map["ID"], id_map["length"]))
                else:
                    print("Error, contig length not match, input vcf may not using b37(hg19) reference. {}".format(line))
                    sys.exit(1)
            elif id_map["ID"] in b37_contigs:
                if int(id_map["length"]) == b37_contigs[id_map["ID"]]:
                    header_new.append(line)
                else:
                    print("Error, contig length not match, input vcf may not using b37(hg19) reference. {}".format(line))
                    sys.exit(1)
        else:
            header_new.append(line)
    return header_new

def rename_sv_sample(header_lines, sample_id):
    header_new = []
    for line in header_lines:
        line = line.strip()
        if line.startswith("#CHROM"):
            fields = line.split("\t")
            fields[9] = sample_id
            header_new.append("\t".join(fields))
        else:
            header_new.append(line)
    return header_new
        
def remove_chr_in_record(record_line):

    record_line = record_line.strip()

    b37_contigs = {'1':249250621, '2':243199373, '3':198022430, '4':191154276, 
        '5':180915260, '6':171115067, '7':159138663, '8':146364022, '9':141213431,
        '10':135534747, '11':135006516, '12':133851895, '13':115169878, '14':107349540,
        '15':102531392, '16':90354753, '17':81195210, '18':78077248, '19':59128983,
        '20':63025520, '21':48129895, '22':51304566, 'X':155270560, 'Y':59373566, 'MT':16569}
    chr_map = dict([("chr"+str(i), str(i)) for i in range(1,23)])
    chr_map["chrX"] = "X"
    chr_map["chrY"] = "Y"
    sv_record = sv_vcf.sv_vcf_record(record_line)
    fields = record_line.split("\t")
    if sv_record.info_dict["SVTYPE"] != "BND":
        if fields[0] in chr_map and sv_record.info_dict["CHR2"] in chr_map:
            fields[0] = chr_map[fields[0]]
            info_fields = fields[7].split(";")
            for i in range(len(info_fields)):
                if info_fields[i].split("=")[0] == "CHR2":
                    info_fields[i] = "CHR2=" + chr_map[sv_record.info_dict["CHR2"]]
            fields[7] = ";".join(info_fields)
            return "\t".join(fields)
        elif fields[0] in b37_contigs and sv_record.info_dict["CHR2"] in b37_contigs:
            return "\t".join(fields)
    else:
        if fields[0] in chr_map and sv_record.chrom2 in chr_map:
            fields[0] = chr_map[fields[0]]
            bnd_string = fields[4]
            fields[4] = bnd_string.replace(sv_record.chrom2, chr_map[sv_record.chrom2])
            if "CHR2=" in fields[7]:
                info_fields = fields[7].split(";")
                for i in range(len(info_fields)):
                    if info_fields[i].split("=")[0] == "CHR2":
                        info_fields[i] = "CHR2=" + chr_map[sv_record.info_dict["CHR2"]]
                fields[7] = ";".join(info_fields)
            return "\t".join(fields)
        elif fields[0] in b37_contigs and sv_record.chrom2 in b37_contigs:
            return "\t".join(fields)
    return ""

def main():
    header_lines = []
    with open(sys.argv[1], "r") as io:
        for line in io:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                break
    sample_id = sys.argv[1].split(".")[0]
    h1 = remove_chr_in_header(header_lines)
    h2 = rename_sv_sample(h1, sample_id)
    
    out_fp = open(sys.argv[2], "w")
    print("\n".join(h2), file=out_fp)
    with open(sys.argv[1], "r") as io:
        for line in io:
            if line.startswith("#"):
                continue
            new_record = remove_chr_in_record(line)
            if new_record != "":
                print(new_record, file=out_fp)

if __name__ == "__main__":
    main()

