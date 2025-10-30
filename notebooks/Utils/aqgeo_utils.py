import numpy as np
from simpeg.directives import SaveEveryIteration,InversionDirective
from simpeg.utils import sdiag,Zero,validate_type
from simpeg.maps import IdentityMap,Wires

def get_ind_prism(hx,hy,hz,phi_x,phi_y,phi_z,x_0,y_0,z_0, xyz):
    """


    """

  # Background conductivity



    # Convert input coordinates toynp arrays
    X = np.array(xyz[:, 0])
    Y = np.array(xyz[:, 1])
    Z = np.array(xyz[:, 2])
    xyz = np.vstack((X, Y, Z))  # Shape: (3, n_points)

    # Process each prism in the model

    # Extract inversion parameters for prism i
    # h are half-widths of the 3D prism, phi are rotation angles (in degrees),
    # x_0, y_0, z_0 are the center coordinates, p_1 is the conductivity value

    # Build the center coordinates as a column vector
    xyz_0 = np.vstack((x_0, y_0, z_0))

    # Create rotation matrix around x-axis
    # Rx rotates coordinates in the y-z plane
    Rx = np.zeros((3, 3))
    Rx[0, 0] = 1
    Rx[1, 1] = np.cos(phi_x * np.pi / 180)  # Convert degrees to radians
    Rx[1, 2] = -np.sin(phi_x * np.pi / 180)
    Rx[2, 2] = Rx[1, 1]  # cos(phix)
    Rx[2, 1] = -Rx[1, 2]  # sin(phix)

    # Create rotation matrix around y-axis
    # Ry rotates coordinates in the x-z plane
    Ry = np.zeros_like(Rx)
    Ry[1, 1] = 1
    Ry[0, 0] = np.cos(phi_y * np.pi / 180)
    Ry[2, 0] = -np.sin(phi_y * np.pi / 180)
    Ry[2, 2] = Ry[0, 0]  # cos(phiy)
    Ry[0, 2] = -Ry[2, 0]  # sin(phiy)

    # Create rotation matrix around z-axis
    # Rz rotates coordinates in the x-y plane
    Rz = np.zeros_like(Rx)
    Rz[2, 2] = 1
    Rz[0, 0] = np.cos(phi_z * np.pi / 180)
    Rz[0, 1] = -np.sin(phi_z * np.pi / 180)
    Rz[1, 1] = Rz[0, 0]  # cos(phiz)
    Rz[1, 0] = -Rz[0, 1]  # sin(phiz)

    # Combine rotation matrices: order matters (Rx @ Ry @ Rz)
    M = Rx @ Ry @ Rz

    # Translate coordinates to prism center and apply rotation
    xyz_m_xyz_0 = xyz - xyz_0  # Shape: (3, n_points)
    tau = (xyz_m_xyz_0).T @ M

    ind_block = (
            (tau[:, 0] > -hx)
            & (tau[:, 0] < hx)
            & (tau[:, 1] > -hy)
            & (tau[:, 1] < hy)
            & (tau[:, 2] > -hz)
            & (tau[:, 2] < hz)
    )

    return ind_block

class UpdatePreconditioner(InversionDirective):
    """
    Create a Jacobi preconditioner for the linear problem
    """

    def __init__(self,solver, update_every_iteration=True,jtj_approx=None,reg_only=False, **kwargs):
        super().__init__(**kwargs)
        self.solver = solver
        self.update_every_iteration = update_every_iteration
        self.jtj_set = jtj_approx
        self.reg_only = reg_only

    @property
    def update_every_iteration(self):
        """Whether to update the preconditioner at every iteration.

        Returns
        -------
        bool
        """
        return self._update_every_iteration

    @update_every_iteration.setter
    def update_every_iteration(self, value):
        self._update_every_iteration = validate_type(
            "update_every_iteration", value, bool
        )

    def initialize(self):
        # Create the pre-conditioner
        regDiag = np.zeros_like(self.invProb.model)
        m = self.invProb.model

        for reg in self.reg.objfcts:
            # Check if regularization has a projection
            rdg = reg.deriv2(m)
            if not isinstance(rdg, Zero):
                regDiag += rdg.diagonal()

        JtJdiag = np.zeros_like(self.invProb.model)
        if not self.reg_only:
            if self.jtj_set is None:
                for sim, dmisfit in zip(self.simulation, self.dmisfit.objfcts):
                    if getattr(sim, "getJtJdiag", None) is None:
                        assert getattr(sim, "getJ", None) is not None, (
                            "Simulation does not have a getJ attribute."
                            + "Cannot form the sensitivity explicitly"
                        )
                        JtJdiag += np.sum(np.power((dmisfit.W * sim.getJ(m)), 2), axis=0)
                    else:
                        JtJdiag += sim.getJtJdiag(m, W=dmisfit.W)
            else:
                JtJdiag = self.jtj_set


        reg =self.invProb.beta * self.reg.deriv2(m)+sdiag(JtJdiag)

        PC=self.solver(reg)


        self.opt.approxHinv = PC

    def endIter(self):
        # Cool the threshold parameter
        if self.update_every_iteration is False:
            return

        # Create the pre-conditioner
        regDiag = np.zeros_like(self.invProb.model)
        m = self.invProb.model

        for reg in self.reg.objfcts:
            # Check if he has wire
            regDiag += reg.deriv2(m).diagonal()

        JtJdiag = np.zeros_like(self.invProb.model)

        if not self.reg_only:
            if self.jtj_set is None:
                for sim, dmisfit in zip(self.simulation, self.dmisfit.objfcts):
                    if getattr(sim, "getJtJdiag", None) is None:
                        assert getattr(sim, "getJ", None) is not None, (
                                "Simulation does not have a getJ attribute."
                                + "Cannot form the sensitivity explicitly"
                        )
                        JtJdiag += np.sum(np.power((dmisfit.W * sim.getJ(m)), 2), axis=0)
                    else:
                        JtJdiag += sim.getJtJdiag(m, W=dmisfit.W)
            else:
                JtJdiag = self.jtj_set

        reg = self.invProb.beta * self.reg.deriv2(m)+sdiag(JtJdiag)

        PC=self.solver(reg)


        self.opt.approxHinv = PC