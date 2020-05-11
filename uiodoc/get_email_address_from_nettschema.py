import pandas as pd
import os
import shutil
from optparse import OptionParser

def main():
    version='v.1.1'
    usage = "usage: %prog --input filename [--output filename]"
    help = "Reads nettschema excel list of registered participants and generates list with email addresses for import in Thunderbird."
    parser = OptionParser(usage=usage)
    parser.add_option("--input", dest="input_file",
                      help="input filename .xls", metavar="input_file",type=str )
    parser.add_option("--output", dest="output_file",
                      help="'output directory'", metavar="path", type=str)
    (options, args) = parser.parse_args()

    if not options.input_file:
        parser.error("An input file must be specified!")
    else:
        input_file = options.input_file
        if input_file.find(".xls") < 0:
            parser.error("Wrong input file type. Input file must be .xls or .xlsx!")
    if not options.output_file:
        output_file = os.path.basename(input_file)[:-5]+".csv"
    else:
        output_file = options.output_file
        
    # Read file
    data = pd.read_excel(input_file)
    columns = ["First Name","Last Name","Display Name","Nickname","Primary Email","Secondary Email","Screen Name","Work Phone","Home Phone","Fax Number","Pager Number","Mobile Number","Home Address","Home Address 2","Home City","Home State","Home ZipCode","Home Country","Work Address","Work Address 2","Work City","Work State","Work ZipCode","Work Country","Job Title","Department","Organization","Web Page 1","Web Page 2","Birth Year","Birth Month","Birth Day","Custom 1","Custom 2","Custom 3","Custom 4","Notes",""]
    # Save data to csv
    csv_out = pd.DataFrame(columns=columns)
    csv_out["Primary Email"] = data["Email address"].values
    csv_out["Display Name"] = data["Email address"].values
    print("Write to file %s" % output_file)
    csv_out.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()
