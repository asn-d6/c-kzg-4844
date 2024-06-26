package ethereum.ckzg4844.test_formats;

import com.fasterxml.jackson.annotation.JsonProperty;
import org.apache.tuweni.bytes.Bytes;

public class VerifyCellKzgProofTest {
  public static class Input {
    private String commitment;

    @JsonProperty("cell_index")
    private Long cellIndex;

    private String cell;
    private String proof;

    public byte[] getCommitment() {
      return Bytes.fromHexString(commitment).toArray();
    }

    public Long getCellIndex() {
      return cellIndex;
    }

    public byte[] getCell() {
      return Bytes.fromHexString(cell).toArrayUnsafe();
    }

    public byte[] getProof() {
      return Bytes.fromHexString(proof).toArray();
    }
  }

  private Input input;
  private Boolean output;

  public Input getInput() {
    return input;
  }

  public Boolean getOutput() {
    return output;
  }
}
