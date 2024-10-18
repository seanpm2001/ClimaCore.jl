import ClimaCore: DataLayouts, Topologies, Spaces, Fields
import ClimaCore.DataLayouts: getindex_field, setindex_field!
using CUDA
import ClimaCore.Topologies
import ClimaCore.Topologies: perimeter_vertex_node_index

_max_threads_cuda() = 256

function _configure_threadblock(max_threads, nitems)
    nthreads = min(max_threads, nitems)
    nblocks = cld(nitems, nthreads)
    return (nthreads, nblocks)
end

_configure_threadblock(nitems) =
    _configure_threadblock(_max_threads_cuda(), nitems)

function Topologies.dss_load_perimeter_data!(
    ::ClimaComms.CUDADevice,
    dss_buffer::Topologies.DSSBuffer,
    data::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    perimeter::Topologies.Perimeter2D,
)
    (; perimeter_data) = dss_buffer
    nitems = prod(size(parent(perimeter_data)))
    args = (perimeter_data, data, perimeter)
    threads = _max_threads_cuda()
    p = linear_partition(nitems, threads)
    auto_launch!(
        dss_load_perimeter_data_kernel!,
        args;
        threads_s = p.threads,
        blocks_s = p.blocks,
    )
    return nothing
end

function dss_load_perimeter_data_kernel!(
    perimeter_data::DataLayouts.AbstractData,
    data::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    perimeter::Topologies.Perimeter2D{Nq},
) where {Nq}
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    (nperimeter, _, _, nlevels, nelems) = size(perimeter_data)
    nfidx = DataLayouts.ncomponents(perimeter_data)
    sizep = (nlevels, nperimeter, nfidx, nelems) # assume VIFH order
    CI = CartesianIndex

    if gidx ≤ prod(sizep)
        (level, p, fidx, elem) = cart_ind(sizep, gidx).I
        (ip, jp) = perimeter[p]
        val = getindex_field(data, CI(ip, jp, fidx, level, elem))
        setindex_field!(perimeter_data, val, CI(p, 1, fidx, level, elem))
    end
    return nothing
end

function Topologies.dss_unload_perimeter_data!(
    ::ClimaComms.CUDADevice,
    data::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    dss_buffer::Topologies.DSSBuffer,
    perimeter,
)
    (; perimeter_data) = dss_buffer
    nitems = prod(size(parent(perimeter_data)))
    args = (data, perimeter_data, perimeter)
    threads = _max_threads_cuda()
    p = linear_partition(nitems, threads)
    auto_launch!(
        dss_unload_perimeter_data_kernel!,
        args;
        threads_s = p.threads,
        blocks_s = p.blocks,
    )
    return nothing
end

function dss_unload_perimeter_data_kernel!(
    data::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    perimeter_data::AbstractData,
    perimeter::Topologies.Perimeter2D{Nq},
) where {Nq}
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    (nperimeter, _, _, nlevels, nelems) = size(perimeter_data)
    nfidx = DataLayouts.ncomponents(perimeter_data)
    sizep = (nlevels, nperimeter, nfidx, nelems) # assume VIFH order
    CI = CartesianIndex

    if gidx ≤ prod(sizep)
        (level, p, fidx, elem) = cart_ind(sizep, gidx).I
        (ip, jp) = perimeter[p]
        val = getindex_field(perimeter_data, CI(p, 1, fidx, level, elem))
        setindex_field!(data, val, CI(ip, jp, fidx, level, elem))
    end
    return nothing
end

function Topologies.dss_local!(
    ::ClimaComms.CUDADevice,
    perimeter_data::DataLayouts.VIFH,
    perimeter::Topologies.Perimeter2D,
    topology::Topologies.Topology2D,
)
    nlocalvertices = length(topology.local_vertex_offset) - 1
    nlocalfaces = length(topology.interior_faces)
    if (nlocalvertices + nlocalfaces) > 0
        (nperimeter, _, _, nlevels, nelems) = size(perimeter_data)
        nfid = DataLayouts.ncomponents(perimeter_data)
        nitems = nlevels * nfid * (nlocalfaces + nlocalvertices)
        args = (
            perimeter_data,
            topology.local_vertices,
            topology.local_vertex_offset,
            topology.interior_faces,
            perimeter,
        )
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)
        auto_launch!(
            dss_local_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
    end
    return nothing
end

function dss_local_kernel!(
    perimeter_data::DataLayouts.VIFH,
    local_vertices::AbstractVector{Tuple{Int, Int}},
    local_vertex_offset::AbstractVector{Int},
    interior_faces::AbstractVector{Tuple{Int, Int, Int, Int, Bool}},
    perimeter::Topologies.Perimeter2D,
)
    FT = eltype(parent(perimeter_data))
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    nlocalvertices = length(local_vertex_offset) - 1
    nlocalfaces = length(interior_faces)
    (nperimeter, _, _, nlevels, _) = size(perimeter_data)
    nfidx = DataLayouts.ncomponents(perimeter_data)
    CI = CartesianIndex
    if gidx ≤ nlevels * nfidx * nlocalvertices # local vertices
        sizev = (nlevels, nfidx, nlocalvertices)
        (level, fidx, vertexid) = cart_ind(sizev, gidx).I
        sum_data = FT(0)
        st, en =
            local_vertex_offset[vertexid], local_vertex_offset[vertexid + 1]
        for idx in st:(en - 1)
            (lidx, vert) = local_vertices[idx]
            ip = perimeter_vertex_node_index(vert)
            sum_data +=
                getindex_field(perimeter_data, CI(ip, 1, fidx, level, lidx))
        end
        for idx in st:(en - 1)
            (lidx, vert) = local_vertices[idx]
            ip = perimeter_vertex_node_index(vert)
            setindex_field!(
                perimeter_data,
                sum_data,
                CI(ip, 1, fidx, level, lidx),
            )
        end
    elseif gidx ≤ nlevels * nfidx * (nlocalvertices + nlocalfaces) # interior faces
        nfacedof = div(nperimeter - 4, 4)
        sizef = (nlevels, nfidx, nlocalfaces)
        (level, fidx, faceid) =
            cart_ind(sizef, gidx - nlevels * nfidx * nlocalvertices).I
        (lidx1, face1, lidx2, face2, reversed) = interior_faces[faceid]
        (first1, inc1, last1) =
            Topologies.perimeter_face_indices_cuda(face1, nfacedof, false)
        (first2, inc2, last2) =
            Topologies.perimeter_face_indices_cuda(face2, nfacedof, reversed)
        for i in 1:nfacedof
            ip1 = inc1 == 1 ? first1 + i - 1 : first1 - i + 1
            ip2 = inc2 == 1 ? first2 + i - 1 : first2 - i + 1
            idx1 = CI(ip1, 1, fidx, level, lidx1)
            idx2 = CI(ip2, 1, fidx, level, lidx2)
            val =
                getindex_field(perimeter_data, idx1) +
                getindex_field(perimeter_data, idx2)
            setindex_field!(perimeter_data, val, idx1)
            setindex_field!(perimeter_data, val, idx2)
        end
    end

    return nothing
end

function Topologies.dss_transform!(
    device::ClimaComms.CUDADevice,
    perimeter_data::DataLayouts.VIFH,
    data::Union{DataLayouts.VIJFH, DataLayouts.IJFH},
    perimeter::Topologies.Perimeter2D,
    local_geometry::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    weight::DataLayouts.IJFH,
    localelems::AbstractVector{Int},
)
    nlocalelems = length(localelems)
    if nlocalelems > 0
        (nperimeter, _, _, nlevels, _) =
            DataLayouts.universal_size(perimeter_data)
        nitems = nlevels * nperimeter * nlocalelems
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)

        args = (
            perimeter_data,
            data,
            perimeter,
            local_geometry,
            weight,
            localelems,
            Val(nlocalelems),
        )
        auto_launch!(
            dss_transform_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
    end
    return nothing
end

function dss_transform_kernel!(
    perimeter_data::DataLayouts.VIFH,
    data::Union{DataLayouts.VIJFH, DataLayouts.IJFH},
    perimeter::Topologies.Perimeter2D,
    local_geometry::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    weight::DataLayouts.IJFH,
    localelems::AbstractVector{Int},
    ::Val{nlocalelems},
) where {nlocalelems}
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    (nperimeter, _, _, nlevels, nelems) =
        DataLayouts.universal_size(perimeter_data)
    CI = CartesianIndex
    if gidx ≤ nlevels * nperimeter * nlocalelems
        sizet = (nlevels, nperimeter, nlocalelems)
        (level, p, localelemno) = cart_ind(sizet, gidx).I
        elem = localelems[localelemno]
        (ip, jp) = perimeter[p]
        loc = CI(ip, jp, 1, level, elem)
        src = Topologies.dss_transform(
            data[loc],
            local_geometry[loc],
            weight[loc],
        )
        perimeter_data[CI(p, 1, 1, level, elem)] =
            Topologies.drop_vert_dim(eltype(perimeter_data), src)
    end
    return nothing
end

function Topologies.dss_untransform!(
    device::ClimaComms.CUDADevice,
    perimeter_data::DataLayouts.VIFH,
    data::Union{DataLayouts.VIJFH, DataLayouts.IJFH},
    local_geometry::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    perimeter::Topologies.Perimeter2D,
    localelems::AbstractVector{Int},
)
    nlocalelems = length(localelems)
    if nlocalelems > 0
        (nperimeter, _, _, nlevels, _) =
            DataLayouts.universal_size(perimeter_data)
        nitems = nlevels * nperimeter * nlocalelems
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)
        args = (
            perimeter_data,
            data,
            local_geometry,
            perimeter,
            localelems,
            Val(nlocalelems),
        )
        auto_launch!(
            dss_untransform_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
    end
    return nothing
end

function dss_untransform_kernel!(
    perimeter_data::DataLayouts.VIFH,
    data::Union{DataLayouts.VIJFH, DataLayouts.IJFH},
    local_geometry::Union{DataLayouts.IJFH, DataLayouts.VIJFH},
    perimeter::Topologies.Perimeter2D,
    localelems::AbstractVector{Int},
    ::Val{nlocalelems},
) where {nlocalelems}
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    (nperimeter, _, _, nlevels, _) = DataLayouts.universal_size(perimeter_data)
    CI = CartesianIndex
    if gidx ≤ nlevels * nperimeter * nlocalelems
        sizet = (nlevels, nperimeter, nlocalelems)
        (level, p, localelemno) = cart_ind(sizet, gidx).I
        elem = localelems[localelemno]
        ip, jp = perimeter[p]

        loc = CI(ip, jp, 1, level, elem)
        data[loc] = Topologies.dss_untransform(
            eltype(data),
            perimeter_data[CI(p, 1, 1, level, elem)],
            local_geometry[loc],
        )
    end
    return nothing
end

# TODO: Function stubs, code to be implemented, needed only for distributed GPU runs
function Topologies.dss_local_ghost!(
    ::ClimaComms.CUDADevice,
    perimeter_data::DataLayouts.VIFH,
    perimeter::Topologies.Perimeter2D,
    topology::Topologies.AbstractTopology,
)
    nghostvertices = length(topology.ghost_vertex_offset) - 1
    if nghostvertices > 0
        (_, _, _, nlevels, _) = size(perimeter_data)
        nfid = DataLayouts.ncomponents(perimeter_data)
        nitems = nlevels * nfid * nghostvertices
        args = (
            perimeter_data,
            topology.ghost_vertices,
            topology.ghost_vertex_offset,
            perimeter,
        )
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)
        auto_launch!(
            dss_local_ghost_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
    end
    return nothing
end

function dss_local_ghost_kernel!(
    perimeter_data::DataLayouts.VIFH,
    ghost_vertices,
    ghost_vertex_offset,
    perimeter::Topologies.Perimeter2D,
)
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    FT = eltype(parent(perimeter_data))
    (nperimeter, _, _, nlevels, _) = size(perimeter_data)
    nfidx = DataLayouts.ncomponents(perimeter_data)
    CI = CartesianIndex
    nghostvertices = length(ghost_vertex_offset) - 1
    if gidx ≤ nlevels * nfidx * nghostvertices
        sizev = (nlevels, nfidx, nghostvertices)
        (level, fidx, vertexid) = cart_ind(sizev, gidx).I
        sum_data = FT(0)
        st, en =
            ghost_vertex_offset[vertexid], ghost_vertex_offset[vertexid + 1]
        for idx in st:(en - 1)
            isghost, lidx, vert = ghost_vertices[idx]
            if !isghost
                ip = perimeter_vertex_node_index(vert)
                sum_data +=
                    getindex_field(perimeter_data, CI(ip, 1, fidx, level, lidx))
            end
        end
        for idx in st:(en - 1)
            isghost, lidx, vert = ghost_vertices[idx]
            if !isghost
                ip = perimeter_vertex_node_index(vert)
                setindex_field!(
                    perimeter_data,
                    sum_data,
                    CI(ip, 1, fidx, level, lidx),
                )
            end
        end
    end
    return nothing
end

function Topologies.fill_send_buffer!(
    ::ClimaComms.CUDADevice,
    dss_buffer::Topologies.DSSBuffer;
    synchronize = true,
)
    (; perimeter_data, send_buf_idx, send_data) = dss_buffer
    (nperimeter, _, _, nlevels, nelems) = size(perimeter_data)
    nfid = DataLayouts.ncomponents(perimeter_data)
    nsend = size(send_buf_idx, 1)
    if nsend > 0
        nitems = nsend * nlevels * nfid
        args = (send_data, send_buf_idx, perimeter_data, Val(nsend))
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)
        auto_launch!(
            fill_send_buffer_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
        if synchronize
            CUDA.synchronize(; blocking = true) # CUDA MPI uses a separate stream. This will synchronize across streams
        end
    end
    return nothing
end

function fill_send_buffer_kernel!(
    send_data::AbstractArray{FT, 1},
    send_buf_idx::AbstractArray{I, 2},
    perimeter_data::AbstractData,
    ::Val{nsend},
) where {FT <: AbstractFloat, I <: Int, nsend}
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    (_, _, _, nlevels, nelems) = size(perimeter_data)
    nfid = DataLayouts.ncomponents(perimeter_data)
    sizet = (nlevels, nfid, nsend)
    CI = CartesianIndex
    if gidx ≤ nlevels * nfid * nsend
        (level, fidx, isend) = cart_ind(sizet, gidx).I
        lidx = send_buf_idx[isend, 1]
        ip = send_buf_idx[isend, 2]
        idx = level + ((fidx - 1) + (isend - 1) * nfid) * nlevels
        send_data[idx] =
            getindex_field(perimeter_data, CI(ip, 1, fidx, level, lidx))
    end
    return nothing
end

function Topologies.load_from_recv_buffer!(
    ::ClimaComms.CUDADevice,
    dss_buffer::Topologies.DSSBuffer,
)
    (; perimeter_data, recv_buf_idx, recv_data) = dss_buffer
    (nperimeter, _, _, nlevels, nelems) = size(perimeter_data)
    nfid = DataLayouts.ncomponents(perimeter_data)
    nrecv = size(recv_buf_idx, 1)
    if nrecv > 0
        nitems = nrecv * nlevels * nfid
        args = (perimeter_data, recv_data, recv_buf_idx, Val(nrecv))
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)
        auto_launch!(
            load_from_recv_buffer_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
    end
    return nothing
end

function load_from_recv_buffer_kernel!(
    perimeter_data::AbstractData,
    recv_data::AbstractArray{FT, 1},
    recv_buf_idx::AbstractArray{I, 2},
    ::Val{nrecv},
) where {FT <: AbstractFloat, I <: Int, nrecv}
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    pperimeter_data = parent(perimeter_data)
    (_, _, _, nlevels, nelems) = size(perimeter_data)
    nfid = DataLayouts.ncomponents(perimeter_data)
    sizet = (nlevels, nfid, nrecv)
    CI = CartesianIndex
    if gidx ≤ nlevels * nfid * nrecv
        (level, fidx, irecv) = cart_ind(sizet, gidx).I
        lidx = recv_buf_idx[irecv, 1]
        ip = recv_buf_idx[irecv, 2]
        idx = level + ((fidx - 1) + (irecv - 1) * nfid) * nlevels
        ci = CI(ip, 1, fidx, level, lidx)
        # CUDA.@atomic has limited support, so
        # let's use the methods in DataLayouts
        # to allow this to work:
        s = DataLayouts.singleton(perimeter_data)
        data_inds = DataLayouts.to_data_specific_field(s, ci.I)
        CUDA.@atomic pperimeter_data[data_inds...] += recv_data[idx]
    end
    return nothing
end


function Topologies.dss_ghost!(
    ::ClimaComms.CUDADevice,
    perimeter_data::DataLayouts.VIFH,
    perimeter::Topologies.Perimeter2D,
    topology::Topologies.Topology2D,
)
    nghostvertices = length(topology.ghost_vertex_offset) - 1
    if nghostvertices > 0
        (_, _, _, nlevels, _) = size(perimeter_data)
        nfidx = DataLayouts.ncomponents(perimeter_data)
        nitems = nlevels * nfidx * nghostvertices
        args = (
            perimeter_data,
            topology.ghost_vertices,
            topology.ghost_vertex_offset,
            topology.repr_ghost_vertex,
            perimeter,
        )
        threads = _max_threads_cuda()
        p = linear_partition(nitems, threads)
        auto_launch!(
            dss_ghost_kernel!,
            args;
            threads_s = p.threads,
            blocks_s = p.blocks,
        )
    end
    return nothing
end

function dss_ghost_kernel!(
    perimeter_data::AbstractData,
    ghost_vertices,
    ghost_vertex_offset,
    repr_ghost_vertex,
    perimeter::Topologies.Perimeter2D,
)
    FT = eltype(parent(perimeter_data))
    gidx = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    (_, _, _, nlevels, _) = size(perimeter_data)
    nfidx = DataLayouts.ncomponents(perimeter_data)
    nghostvertices = length(ghost_vertex_offset) - 1
    CI = CartesianIndex
    if gidx ≤ nlevels * nfidx * nghostvertices
        (level, fidx, ghostvertexidx) =
            cart_ind((nlevels, nfidx, nghostvertices), gidx).I
        idxresult, lvertresult = repr_ghost_vertex[ghostvertexidx]
        ipresult = perimeter_vertex_node_index(lvertresult)
        result = getindex_field(
            perimeter_data,
            CI(ipresult, 1, fidx, level, idxresult),
        )
        st, en = ghost_vertex_offset[ghostvertexidx],
        ghost_vertex_offset[ghostvertexidx + 1]
        for vertexidx in st:(en - 1)
            isghost, eidx, lvert = ghost_vertices[vertexidx]
            if !isghost
                ip = perimeter_vertex_node_index(lvert)
                setindex_field!(
                    perimeter_data,
                    result,
                    CI(ip, 1, fidx, level, eidx),
                )
            end
        end
    end
    return nothing
end
