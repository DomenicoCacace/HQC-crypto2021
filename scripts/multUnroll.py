import sys
MASKS = int(sys.argv[1])

def shares_size(i, qty):
    if(i < MASKS-1):
        return qty+"/"+str(MASKS)
    return qty + " - (" + qty + "/" + str(MASKS) + ")*" + str(i)

print("// MULTIPLICATION - PART 1")
for i in range(0, MASKS):
    print("memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);")
    print("fast_convolution_mult(raw_temp+("+str(i)+"*(VEC_N_SIZE_64/" + str(MASKS) + ")),\n\t\t\t\t\t  a1+("+str(i)+"*(weight/" + str(MASKS) +")), a2+("+str(i)+"*(VEC_N_SIZE_64/" + str(MASKS) + ")),\n\t\t\t\t\t  "+ shares_size(i, "weight") + ", " + shares_size(i, "VEC_N_SIZE_64") + ");")
    print("reduce(o->s"+str(i)+", raw_temp);")

print("// MULTIPLICATION - PART 2")
for i in range(0, MASKS):
    for j in range(i+1, MASKS):
        print("shake_prng(seed, SEED_BYTES);")
        print("seedexpander_init(&mask_seedexpander, seed, SEED_BYTES);")
        print("vect_set_random_fixed_weight(&mask_seedexpander, s, weight);")
        print("memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);")
        print("fast_convolution_mult(raw_temp+("+str(j)+"*(VEC_N_SIZE_64/" + str(MASKS) + ")),\n\t\t\t\t\t  a1+("+str(i)+"*(weight/" + str(MASKS) +")), a2+("+str(j)+"*(VEC_N_SIZE_64/" + str(MASKS) + ")),\n\t\t\t\t\t  " + shares_size(i, "weight") + ", " + shares_size(j, "VEC_N_SIZE_64") + ");")
        print("reduce(temp1, raw_temp);")
        print("vect_add(temp1, temp1, s, VEC_N_SIZE_64);")
        print("memset(raw_temp, 0x00, (VEC_N_SIZE_64*2+1)*8);")
        print("fast_convolution_mult(raw_temp+("+str(i)+"*(VEC_N_SIZE_64/" + str(MASKS) + ")),\n\t\t\t\t\t  a1+("+str(j)+"*(weight/" + str(MASKS) +")), a2+("+str(i)+"*(VEC_N_SIZE_64/" + str(MASKS) + ")),\n\t\t\t\t\t  " + shares_size(j, "weight") + ", " + shares_size(i, "VEC_N_SIZE_64") + ");")
        print("reduce(temp2, raw_temp);")
        print("vect_add(s1, temp1, temp2, VEC_N_SIZE_64);")
        print("vect_add(o->s"+str(i)+", o->s"+str(i)+", s, VEC_N_SIZE_64);")
        print("vect_add(o->s"+str(j)+", o->s"+str(j)+", s1, VEC_N_SIZE_64);")
        print()

