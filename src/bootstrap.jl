#! /usr/bin/julia

###############################################################################
###############   bootrstrap method to obtain new normal modes    #############
#   method by Sandra Maguid & Sebastián Fernandez Alberti
#   code by pgbarletta
#########################################################
# Safety pig included:
#
#    _._ _..._ .-',     _.._(`))
#   '-. `     '  /-._.-'    ',/
#      )         \            '.
#     / _    _    |             \
#    |  a    a    /              |
#    \   .-.                     ;  
#     '-('' ).-'       ,'       ;
#        '-;           |      .'
#           \           \    /
#           | 7  .__  _.-\   \
#           | |  |  ``/  /`  /
#          /,_|  |   /,_/   /
#             /,_/      '`-'
###############################################################################
using StatsBase
using DataFrames
using Formatting
##########
# functions
##########
function tognm(vtor_anm)
    vtor_gnm = Array{Float64}(convert(Int64, length(vtor_anm)/3));
    vtor_anm =  vtor_anm.^2
    for i=1:convert(Int64, length(vtor_anm)/3)
        vtor_gnm[i] = sqrt(vtor_anm[i*3-2] + vtor_anm[i*3-1] + vtor_anm[i*3])
    end
    return vtor_gnm
end
#########
function vtor_correl(vtor_input)
    # Preparo variables
    vtor_in = copy(vtor_input)
    vtor_in_corr = Array{Float64}(1);
    vtor_in_corr[1] = 5.0
    marca = length(vtor_in)
    
    for i=1:marca
        # Saco 2 vectores del vtor_input. A uno le saco el 1er elemento y a otro el último        
        vtor_in_prev = copy(vtor_input[1:end-i])
        shift!(vtor_in)
        # Obtengo correlacion entre los 2 vtores (autocorrelacion del vtor original) y determino
        # si es suficientemente alta como p/ seguir calculando correlacion
        corr = corspearman(vtor_in_prev, vtor_in)
        push!(vtor_in_corr, corr)
        if corr < exp(-1)
            marca = copy(i)
            break
        end
    end
    shift!(vtor_in_corr)
    # Devuelvo la secuencia de correlaciones q tuve y el nro de elementos q pude sacar antes de
    # q la correlación baje demasiado
    return vtor_in_corr, marca
end
##########
function writeFortranMatrix(mtx_filename, mtx_to_write)
    mtx_io = open(mtx_filename, "w")
    FORTRAN_format = FormatExpr("{1:< .4E} ")
    FORTRAN_newline_format = FormatExpr("\n")

    for i=1:size(mtx_to_write)[1]
        for j=1:size(mtx_to_write)[2]
            printfmt(mtx_io, FORTRAN_format, mtx_to_write[i,j])
        end    
        printfmt(mtx_io, FORTRAN_newline_format, "")
    end
end

##########
# main program 
##########
# Tomo argumentos de consola
message = string("\n\nUsage:\n",  "./bootstrap.jl <input matrix> ",
"<output matrices prefix> <number of matrices to generate> \"name of info file\"",
"\n\n")
if length(ARGS) < 3 || length(ARGS) > 4
    throw(ArgumentError(message))
elseif length(ARGS) == 4
    output_info_filename = ARGS[4]
    global output_info_bool = true
else
    output_info_bool = false
end
modos_filename = ARGS[1]
out_boot_mtx_filename = ARGS[2]
boot_mtx_num = ARGS[3]
boot_mtx_num = parse(Int64, boot_mtx_num)

# Leo la mtx de modos a bootstrappear
modos = readdlm(modos_filename)
# Preparo variables
corr_mtx = zeros(modos)
autocorr = Array{Int64}(size(modos)[2])
boot_mtx = Array{Float64}(size(modos))

# Obtengo la autocorrelación de los vectores de los modos
for i=1:size(modos)[2]
# Paso el vector ANM a GNM
    gvec = tognm(modos[:, i])
    vec_corr, autocorr[i] = vtor_correl(gvec)
    corr_mtx[1:length(vec_corr), i] = vec_corr
    
end
# Recorto la matriz con las autocorrelaciones 
corr_mtx = corr_mtx[1:maximum(autocorr), :]
# El tamaño de los bloques está hecho p/ GNM, lo adaptop p/ ANM
block = autocorr .* 3;
# Número de bloques del tamaño correspondiente. C/ vector tiene otro bloque extra, ya q rara vez
# son divisibles por el tamaño de bloque q le corresponde.
n_block = convert(Array{Int64}, floor(size(modos)[1] ./ block))
n_block = n_block .+ 1
# Acá obtengo el resto q no queda cubierto por los bloques
resto = size(modos)[1] .% block;

##### Preparo variables y generelo las boot_mtx_num pedidas
rn = Array{Float64}(maximum(n_block))

for mat=1:boot_mtx_num
    out_boot_mtx_filename_current = string(out_boot_mtx_filename, "_", mat)

    for i=1:size(modos)[2]
    # Array con los distintos bloques de este vector q voy a ir eligiendo aleatoriamente
            sample = collect(1:n_block[i])
    # Inicializo los indices inferior y superior de los bloques de la nueva mtx bootstrappeada
            boot_lower = 1
            boot_upper = boot_lower + block[i] - 1
    
            for j=1:n_block[i]
    # Obtengo aleatoriamente el nro de bloque de la columna i de la mtx ingresada q voy a poner en la columna i
    # de la mtx de bootstrap
                rn[j] = rand(sample)
                if rn[j] == n_block[i]
    # Si es el último bloque, tengo q hacer algo particular
                    orig_upper = size(modos)[1]
                    orig_lower = orig_upper - resto[i] + 1
                    boot_upper = boot_lower + resto[i] - 1
    
                    boot_mtx[boot_lower:boot_upper, i] = modos[orig_lower:orig_upper, i]
        
                    boot_lower = boot_lower + resto[i]
                    boot_upper = boot_lower + block[i] - 1
    # Saco el bloque recién copiado de el array 'sample' p/ no volver a copiarlo
                    usado = find(x -> x==rn[j], sample)
                    deleteat!(sample, usado)
                    continue
                end
    # Determino q bloque de la original se copia, lo copio y avanzo al proximo bloque de la bootstrap    
                orig_upper = convert(Int64, rn[j] * block[i])
                orig_lower = convert(Int64, orig_upper - block[i] + 1)
         
                boot_mtx[boot_lower:boot_upper, i] = modos[orig_lower:orig_upper, i]
        
                boot_lower = boot_lower + block[i]
                boot_upper = boot_lower + block[i] - 1
    
    # Saco el bloque recién copiado de el array 'sample' p/ no volver a copiarlo
                usado = find(x -> x==rn[j], sample)
                deleteat!(sample, usado);
    
        end
    end
        
    # Escribo la mtx bootstrappeada
    writeFortranMatrix(out_boot_mtx_filename_current, boot_mtx)    
end
if output_info_bool == true
    writeFortranMatrix(output_info_filename, corr_mtx)
end
