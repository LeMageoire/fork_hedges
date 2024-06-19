# HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints
# William H. Press, John A. Hawkins, Stephen Knox Jones Jr, Jeffrey M. Schaub, Ilya J. Finkelstein
# submitted to Proceedings of the National Academy of Sciences
#
# Demonstration driver and installation validation program
#
# This file demonstrates the use of the separately compiled modules NRpyDNAcode and NRpyRS.
# (Those modules are compiled in C++ using Numerical Recipes.  See separate documentation.)
#
# We encode a specified number of packets from known plaintext, create DNA errors, then
# decode the DNA and compare.

from numpy import *
import NRpyDNAcode as code
import NRpyRS as RS
import argparse

# functions to create sequential packets from the plaintext source, and R-S protect them
def createmesspacket(packno, strandsperpacket, bytesperstrand, messbytesperstrand, strandsperpacketmessage, strandIDbytes, getwiz, wiz_state) : # packno in range 0..255 with value 2 for strandIDbytes
    packet = zeros([strandsperpacket,bytesperstrand],dtype=uint8)
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand,dtype=uint8)
    eof_reached = False
    for i in range(strandsperpacket) :
        packet[i,0] = packno # note assumes value 2 for strandIDbytes
        packet[i,1] = i
        if i < strandsperpacketmessage and not eof_reached :
            ptext = getwiz(messbytesperstrand, wiz_state)
            packet[i,strandIDbytes:strandIDbytes+messbytesperstrand] = ptext
            plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = ptext
            if wiz_state['offset'] == wiz_state['length']:
                eof_reached = True  # Stop filling further strands with data
    return (packet,plaintext)

def protectmesspacket(packetin, strandsperpacket, messbytesperstrand, strandIDbytes) : # fills in the RS check strands
    
    packet = packetin.copy()
    regin = zeros(strandsperpacket,dtype=uint8)
    for j in range(messbytesperstrand) :
        for i in range(strandsperpacket) :
            regin[i] = packet[i,((j+i)% messbytesperstrand)+strandIDbytes]
        regout = RS.rsencode(regin)
        for i in range(strandsperpacket) :
            packet[i,((j+i)% messbytesperstrand)+strandIDbytes] = regout[i]
    return packet

# functions to encode a packet to DNA strands, and decode DNA strands to a packet
def messtodna(mpacket, strandsperpacket, totstrandlen, rightlen) :
    # HEDGES encode a message packet into strands of DNA
    filler = array([0,2,1,3,0,3,2,1,2,0,3,1,3,1,2,0,2,3,1,0,3,2,1,0,1,3],dtype=uint8)
    dpacket = zeros([strandsperpacket,totstrandlen],dtype=uint8)
    for i in range(strandsperpacket) :
        dna = code.encode(mpacket[i,:])
        if len(dna) < totstrandlen : # need filler after message and before right primer
            dnaleft = dna[:-rightlen]
            dnaright = dna[-rightlen:]
            dna = concatenate((dnaleft,filler[:totstrandlen-len(dna)],dnaright))
            #n.b. this can violate the output constraints (very slightly at end of strand)
        elif len(dna) > totstrandlen : 
            raise ValueError("DNA strand too long")
        dpacket[i,:len(dna)] = dna
    return dpacket


def dnatomess(dnapacket, strandsperpacket, bytesperstrand, messbytesperstrand) :
    # HEDGES decode strands of DNA (assumed ordered by packet and ID number) to a packet
    baddecodes = 0
    erasures = 0
    mpacket = zeros([strandsperpacket,bytesperstrand],dtype=uint8)
    epacket = ones([strandsperpacket,bytesperstrand],dtype=uint8) # everything starts as an erasure
    for i in range(strandsperpacket) :
        (errcode, mess, _, _, _, _) = code.decode(dnapacket[i,:],8*bytesperstrand)
        if errcode > 0 :
            baddecodes += 1
            erasures += max(0,messbytesperstrand-len(mess))
        lenmin = min(len(mess),bytesperstrand)
        mpacket[i,:lenmin] = mess[:lenmin]
        epacket[i,:lenmin] = 0
    return (mpacket,epacket,baddecodes,erasures)

#functions to R-S correct a packet and extract its payload to an array of bytes
def correctmesspacket(packetin,epacket, strandsperpacket, messbytesperstrand, strandIDbytes) :
    # error correction of the outer RS code from a HEDGES decoded packet and erasure mask
    packet = packetin.copy()
    regin = zeros(strandsperpacket,dtype=uint8)
    erase = zeros(strandsperpacket,dtype=uint8)
    tot_detect = 0
    tot_uncorrect = 0
    max_detect = 0
    max_uncorrect = 0
    toterrcodes = 0
    for j in range(messbytesperstrand) :
        for i in range(strandsperpacket) :
            regin[i] = packet[i,((j+i)% messbytesperstrand)+strandIDbytes]
            erase[i] = epacket[i,((j+i)% messbytesperstrand)+strandIDbytes]
        locations = array(argwhere(erase),dtype=int32)
        (decoded, errs_detected, errs_corrected, errcode, ok) = RS.rsdecode(regin,locations)
        tot_detect += errs_detected
        tot_uncorrect += max(0,(errs_detected-errs_corrected))
        max_detect = max(max_detect,errs_detected)
        max_uncorrect = max(max_uncorrect,max(0,(errs_detected-errs_corrected)))
        toterrcodes += (0 if errcode==0 else 1)
        for i in range(strandsperpacket) :
            packet[i,((j+i)% messbytesperstrand)+strandIDbytes] = decoded[i]
    return (packet,tot_detect,tot_uncorrect,max_detect,max_uncorrect,toterrcodes)


def extractplaintext(cpacket, strandsperpacketmessage, messbytesperstrand, strandIDbytes) :
    # extract plaintext from a corrected packet
    plaintext = zeros(strandsperpacketmessage*messbytesperstrand,dtype=uint8)
    for i in range(strandsperpacketmessage) :
        plaintext[i*messbytesperstrand:(i+1)*messbytesperstrand] = (
            cpacket[i,strandIDbytes:strandIDbytes+messbytesperstrand])
    return plaintext

# function to create errors in a bag (or packet) of DNA strands
def createerrors(dnabag,srate,drate,irate) :
    # for testing: create errors in a bag of strands
    (nrows,ncols) = dnabag.shape
    newbag = zeros([nrows,ncols],dtype=uint8)
    for i in range(nrows) :
        dna = code.createerrors(dnabag[i,:],srate,drate,irate)
        lenmin = min(len(dna),ncols)
        newbag[i,:lenmin] = dna[:lenmin]
    return newbag

# for every strand in the packet, print the strand data for debugging
def print_strand_data_int(packet, strandsperpacket, messbytesperstrand, strandIDbytes):
    print ("Printing strand data for debugging:")
    for i in range(strandsperpacket):
        strand_data = packet[i, strandIDbytes:strandIDbytes + messbytesperstrand]
        print("Strand {}: {}".format(i, list(strand_data))) # Convert to list for easier reading

def print_strand_data_char(packet, strandsperpacket, messbytesperstrand, strandIDbytes):
    print ("Printing strand data as strings for debugging:")
    for i in range(strandsperpacket):
        strand_data = packet[i, strandIDbytes:strandIDbytes + messbytesperstrand]
        strand_string = ''.join(chr(num) for num in strand_data)  # Convert each byte to character
        print("Strand {}: {}".format(i, strand_string))


def print_all(strandsperpacketcheck, strandsperpacket, file_size_bytes, totstrandlen, strandlen, coderatecode, coderates, bytesperstrand, messbytesperpacket, npackets, messbytesperstrand):
    print("-------------------")
    print("File size: {} bytes".format(file_size_bytes))  # bytes of "D" file
    print("Total strand length: {}".format(totstrandlen))
    print("strandlen: {}".format(strandlen))              # 300 defined in the paper ? (given in original code)
    print("Coderate: {}".format(coderates[coderatecode])) # 0.5 => 50%
    print("Bytes per strand: {}".format(bytesperstrand))  #it's 31.75 => 31 bytes
    print("-------------------")
    print("message bytes per strand: {}".format(messbytesperstrand)) # removed 2 bytes for ID and 2 bytes for runout => 27 bytes
    print("strand per packet: {}".format(255)) # 255 strands per packet (given in original code because of RS(255,32) code)
    print("strand per packet check: {}".format(32)) # 32 strands per packet check (given in original code because of RS(255,32) code)
    print("strand per packet message: {}".format(255-32)) # 223 strands per packet message (given in original code because of RS(255,32) code)
    print("-------------------")
    print("Payload bytes per packet: {}".format(messbytesperpacket))
    payload_bits_per_packet = messbytesperstrand * 8 * (strandsperpacket - strandsperpacketcheck)
    total_nucleotides_per_packet = strandsperpacket * totstrandlen
    bitrate = float(payload_bits_per_packet) / total_nucleotides_per_packet
    print("Packet bitrate: {:.2f} bit/nt".format(bitrate))
    print("Packet ntrate: {:.2f} nt/bit".format(1.0 / bitrate))
    print("Number of packets needed: {}".format(npackets)) # it should be 1 packet for 4878 bytes < X bytes
    print("-------------------")
    print("-------------------")


def main():   
    """
    coderatecode is used to select the code rate from the table of coderates (for ?)
    npackets is the number of packets to generate and test (255 strands each)
    totstrandlen is the total length of DNA strand (300 ?)
    """

    parser = argparse.ArgumentParser(description="Provide substitution, deletion and insertion error rates to test the HEDGES error-correcting code for DNA storage.")
    parser.add_argument('--sub', '-s', dest='e_sub', type=float, action='store', help="substitution error rate", required=True)
    parser.add_argument('--del', '-d', dest='e_del', type=float, action='store', help="deletion error rate", required=True)
    parser.add_argument('--ins', '-i', dest='e_ins', type=float, action='store', help="insertion error rate", required=True)
    args = parser.parse_args()
    # Set parameters
    # DO NOT CHANGE THESE PARAMETERS because it will break the code
    coderates = array([NaN, 0.75, 0.6, 0.5, 1./3., 0.25, 1./6.])  # table of coderates 1..6
    coderatecode = 3  # test this coderate in coderates table above
    #npackets = 20  # number of packets (of 255 strands each) to generate and test
    totstrandlen = 300  # total length of DNA strand
    strandIDbytes = 2  # ID bytes each strand for packet and sequence number
    strandrunoutbytes = 2  # confirming bytes end of each strand (see paper)
    hlimit = 1000000  # maximum size of decode heap, see paper
    leftprimer = "TCGAAGTCAGCGTGTATTGTATG" # 23 nt left primer
    rightprimer = "TAGTGAGTGCGATTAAGCGTGTT"# 23 nt right primer for direct right appending (no revcomp)


    file_size_bytes = 4878
    #file_size_bytes = 15045
    #npackets = ceil(file_size_bytes / messbytesperpacket)

    # Error rates
    (srate, drate, irate) = 1.0 * array([args.e_sub, args.e_del, args.e_ins])
    #print("Substitution error rate: {}".format(srate))
    #print("Deletion error rate: {}".format(drate))
    #print("Insertion error rate: {}".format(irate))
    # DNA constraints
    max_hpoly_run = 4  # max homopolymer length allowed (0 for no constraint)
    GC_window = 12  # window for GC count (0 for no constraint)
    max_GC = 8  # max GC allowed in window (0 for no constraint)
    min_GC = GC_window - max_GC

    # Set additional derived parameters
    leftlen = len(leftprimer)
    rightlen = len(rightprimer)
    strandlen = totstrandlen - leftlen - rightlen
    strandsperpacket = 255  # for RS(255,32)
    strandsperpacketcheck = 32  # for RS(255,32)
    strandsperpacketmessage = strandsperpacket - strandsperpacketcheck

    # Initialize parameters in NRpyDNAcode module
    (NSALT, MAXSEQ, NSTAK, HLIMIT) = code.getparams()  # get settable code parameters
    code.setparams(8*strandIDbytes, MAXSEQ, NSTAK, hlimit)  # change NSALT and HLIMIT
    bytesperstrand = int(strandlen * coderates[coderatecode] / 4.)
    messbytesperstrand = bytesperstrand - strandIDbytes - strandrunoutbytes  # payload bytes per strand
    messbytesperpacket = strandsperpacketmessage * messbytesperstrand  # payload bytes per packet of 255 strands
    code.setcoderate(coderatecode, leftprimer, rightprimer)  # set code rate with left and right primers
    code.setdnaconstraints(GC_window, max_GC, min_GC, max_hpoly_run)  # set DNA constraints (see paper)

    npackets = int(ceil(float(file_size_bytes) / messbytesperpacket))

   # Using .format() for string formatting
    print_all(strandsperpacketcheck,
              strandsperpacket,
              file_size_bytes,
              totstrandlen,
              strandlen, 
              coderatecode, 
              coderates, 
              bytesperstrand, 
              messbytesperpacket, 
              npackets, 
              messbytesperstrand)

    # Define a source of plaintext bytes
    UseWiz = True
    wiz_state = {'offset': 0, 'length': 0, 'bytes': None}
    if UseWiz:
        wizfile = "data/D"
        with open(wizfile, 'r') as myfile:
            wiztext = myfile.read()
        wiz_state['bytes'] = array([c for c in wiztext]).view(uint8)
        wiz_state['length'] = len(wiz_state['bytes'])
        def getwiz(n, wiz_state):
            if wiz_state['offset'] + n > wiz_state['length']:
                remaining = wiz_state['length'] - wiz_state['offset']
                result = wiz_state['bytes'][wiz_state['offset']:wiz_state['length']]
                # Fill the rest with zeros (or another stop code) to denote no more data
                result = array(list(result) + [0] * (n - remaining), dtype=uint8)
                wiz_state['offset'] = wiz_state['length']  # Move offset to end to stop further reads
            else:
                result = wiz_state['bytes'][wiz_state['offset']:wiz_state['offset'] + n]
                wiz_state['offset'] += n
            return result
    else:
        def getwiz(n):
            return random.randint(0, high=256, size=n, dtype=uint8)

    # Script starts executing here
    print("")
    print("Starting encoding and decoding test...")
    print ("HEDGES Error-Correcting Code for DNA Storage Corrects Indels and Allows Sequence Constraints")
    print ("for each packet, these statistics are shown in two groups:")
    print ("1.1 HEDGES decode failures, 1.2 HEDGES bytes thus declared as erasures")
    print ("1.3 R-S total errors detected in packet, 1.4 max errors detected in a single decode")
    print ("2.1 R-S reported as initially-uncorrected-but-recoverable total, 2.2 same, but max in single decode")
    print ("2.3 R-S total error codes; if zero, then R-S corrected all errors")
    print ("2.4 Actual number of byte errors when compared to known plaintext input")
    print ("")


    badpackets = 0
    Totalbads = zeros(8,dtype=int)
    for ipacket in range(npackets) :
        # encode
        messpack, messplain =  createmesspacket(ipacket, strandsperpacket, bytesperstrand, messbytesperstrand, strandsperpacketmessage, strandIDbytes, getwiz, wiz_state) # plaintext to message packet
        #print_strand_data_char(messpack, strandsperpacket, messbytesperstrand, strandIDbytes)
        rspack = protectmesspacket(messpack, strandsperpacket, messbytesperstrand, strandIDbytes) # Reed-Solomon protect the packet
        dnapack = zeros([strandsperpacket,totstrandlen],dtype=uint8)
        try :
            dnapack = messtodna(rspack, strandsperpacket, totstrandlen, rightlen) # encode to strands of DNA containing payload messplain
        except (ValueError) as e :
            print("Error: DNA strand too long -", str(e))
            exit(1)
        # simulate errors in DNA synthesis and sequencing
        obspack = createerrors(dnapack,srate,drate,irate)

        # decode
        (dpacket,epacket,baddecodes,erasures) = dnatomess(obspack, strandsperpacket, bytesperstrand, messbytesperstrand) # decode the strands
        (cpacket,tot_detect,tot_uncorrect,max_detect,max_uncorrect,toterrcodes) = correctmesspacket(dpacket,epacket, strandsperpacket, messbytesperstrand, strandIDbytes)
    
        # check against ground truth
        messcheck = extractplaintext(cpacket, strandsperpacketmessage, messbytesperstrand, strandIDbytes)
        with open("results/results.txt", "wb") as file:
            file.write(messcheck)
        badbytes = count_nonzero(messplain-messcheck)
        # print summary line
        Totalbads += array([ baddecodes,erasures,tot_detect,max_detect,tot_uncorrect,max_uncorrect,toterrcodes,badbytes])
        print ("%3d: (%3d %3d %3d %3d) (%3d %3d %3d %3d)" % (ipacket, baddecodes,erasures,
            tot_detect,max_detect, tot_uncorrect,max_uncorrect,toterrcodes,badbytes)),
        print ("packet OK" if badbytes == 0 else "packet NOT ok")
        if badbytes : badpackets += 1
    #print ("all packets OK" if not badpackets else "some packets had errors!")
    #print ("TOT: (%4d %4d %4d %4d) (%4d %4d %4d %4d)" % tuple(Totalbads))
    if badpackets :
        print ("Some packets had errors")
        exit(1)
    else :
        print ("All packets OK")
        exit(0)

if __name__ == "__main__":
    main()