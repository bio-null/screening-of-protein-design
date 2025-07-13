import argparse


'''
该文件接受在命令行接受用户输入的参数，以改变proteinmpnn设计蛋白质产生的fa文件中的新序列的命名方式，以方便

'''




def parsed_args():
    parser=argparse.ArgumentParser(description="rename the sequences in fasta file")
    parser.add_argument('-i','--input',required=True,help="the fasta file to rename")
    parser.add_argument('-o','--output',help="file for data to be dumped,default=inputed file")
    parser.add_argument('-p','--prefix',help="prefix to rename the sequence,default=protein 's name")
    return parser.parse_args()


def main():
    args=parsed_args()
    input_path=args.input
    output_path=args.output if args.output else args.input
    with open(input_path,'r') as file1,open(output_path,'w') as file2:
        line1=file1.readline()
        lin2=file1.readline()
        protein_name=line1.split(',')[0][1:]
        prefix=args.prefix if args.prefix else protein_name
        for line in file1:
            if line.startswith(">"):
                words=line.split(",")
                words[0]=">"+prefix+words[1].split("=")[-1]
                for word in words[:-1]:
                    file2.write(f"{word},")
                file2.write(f"{words[-1]}")
                file2.write("\n")
            else:
                file2.write(f"{line}\n")


if __name__=="__main__":
    main()
